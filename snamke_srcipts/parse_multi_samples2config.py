import re
import pandas as pd
import yaml
import argparse

from pathlib import Path
from collections import OrderedDict, defaultdict

from yaml.representer import Representer 
yaml.add_representer(defaultdict, Representer.represent_dict)
yaml.add_representer(OrderedDict, Representer.represent_dict)

def read_table(file):
    suffix = Path(file).suffix
    assert suffix in ['.txt', '.xlsx', '.csv'], f'suffix: {suffix},Input file should be .txt/.csv/.xlsx format.'
    if suffix == '.txt':
        return pd.read_table(file)
    elif suffix == '.xlsx':
        return pd.read_excel(file)
    else:
        return pd.read_csv(file)
    
def reverse_complement(seq):
    return seq.upper().translate(str.maketrans('ATCG', 'TAGC'))[::-1]

def get_tag_seq(primer, insert_tags, uplib=True, kmer: int=2):
    if not uplib:
        insert_tags = reverse_complement(insert_tags)
        # insert_tags = insert_tags[::-1]
    try:
        extend = re.search(rf'(^.*{primer})(\w\w)', insert_tags).group(2)
    except:
        raise KeyError(f'primer seq {primer} and insert tag seq {insert_tags} is not paired, please check.')
    return extend

def extract_adapter_seq(insert_seq: str, tn5_adapter: str, adpater_names: tuple, kmer=15) -> dict: 
    # kmer should be suitable for your experiment
    adapters = {adpater_names[0]: {}, adpater_names[1]: {}}
    insert_seq_rc = reverse_complement(insert_seq)
    tn5_adapter_rc = reverse_complement(tn5_adapter)
    #
    plus, minus = {}, {}
    plus['five_prime_end_read1_adapter'] = tn5_adapter
    plus['five_prime_end_read2_adapter'] = insert_seq[:kmer]
    plus['three_prime_end_read1_adapter'] = insert_seq_rc[-kmer:]
    plus['three_prime_end_read2_adapter'] = tn5_adapter_rc
    minus['five_prime_end_read1_adapter'] = tn5_adapter
    minus['five_prime_end_read2_adapter'] = insert_seq_rc[-kmer-3-3:-3]
    minus['three_prime_end_read1_adapter'] = insert_seq[3:3+kmer+3]
    minus['three_prime_end_read2_adapter'] = tn5_adapter_rc
    adapters[adpater_names[0]].update(plus)
    adapters[adpater_names[1]].update(minus)
    return adapters

def generate_config(df, args):
    config = OrderedDict()

    config['Genome'] = "GRCh38" #基因组序列文件前缀,eg: hg38
    config['Genomedir'] = "/path/to/hg38" # 基因序列文件路径
    config['ChromSize'] = "/path/to/hg38/GRCh38.chrom.sizes"  # 基因组每条染色体长度信息,格式 = 染色体名称<tab>长度

    config['Bin'] = "/path/to/PEGUIDE" # peguide软件路径

    assert len(df['F_lib'].unique()) == len(df['R_lib'].unique()) == 1, f"Noticing the lib name must be unique: F_lib: {df['F_lib'].unique()} R_lib: {df['R_lib'].unique()}"
    config['primer_F'] = df['F_lib'][0] # 正向测序名称
    config['primer_R'] = df['R_lib'][0] # 反向测序名称

    # primer extension sequence
    # config['primer_F_seq'] = "TC"
    # config['primer_R_seq'] = "GG"

    # 软件IGV安装路径
    igv_path = "/path/to/igv"
    config['igvtools'] = f"java -Xmx1500m --module-path={igv_path}/lib @{igv_path}/igv.args --module=org.igv/org.broad.igv.tools.IgvTools"

    config['thread'] = args.thread

    config['Samples'] = list(df['sample_name'].values)

    config['Targets'] = defaultdict(dict)
    config['adapters'] = defaultdict(dict)
    config['extend_seq'] = defaultdict(dict)
    for i in range(df.shape[0]):
        sample = df.loc[i, 'sample_name']
        target_name = df.loc[i, 'target_name'] if 'target_name' in df.columns else sample
        target_seq = df.loc[i, 'target_seq']
        config['Targets'][sample]['name'] = target_name
        config['Targets'][sample]['seq'] = target_seq
        # adapters info
        # config['adapters'][sample] = defaultdict(dict)
        insert_tag = df.loc[i, 'insert_tag']
        combine_F_name = f"{sample}.{config['primer_F']}"
        combine_R_name = f"{sample}.{config['primer_R']}"
        if 'Tn5_adaptor' in df.columns:
            tn5_adapter = df.loc[i, 'Tn5_adaptor']
            info = extract_adapter_seq(insert_tag, tn5_adapter, (combine_F_name, combine_R_name), kmer=args.kmer)
            config['adapters'].update(info)
        else:
            config['adapters'][combine_F_name]['five_prime_end_read1_adapter'] = df.loc[i, 'F_lib_R1_5prime_adaptor_seq']
            config['adapters'][combine_F_name]['five_prime_end_read2_adapter'] = df.loc[i, 'F_lib_R2_5prime_adaptor_seq']
            config['adapters'][combine_F_name]['three_prime_end_read1_adapter'] = df.loc[i, 'F_lib_R1_3prime_adaptor_seq']
            config['adapters'][combine_F_name]['three_prime_end_read2_adapter'] = df.loc[i, 'F_lib_R2_3prime_adaptor_seq']
            # 
            config['adapters'][combine_R_name]['five_prime_end_read1_adapter'] = df.loc[i, 'R_lib_R1_5prime_adaptor_seq']
            config['adapters'][combine_R_name]['five_prime_end_read2_adapter'] = df.loc[i, 'R_lib_R2_5prime_adaptor_seq']
            config['adapters'][combine_R_name]['three_prime_end_read1_adapter'] = df.loc[i, 'R_lib_R1_3prime_adaptor_seq']
            config['adapters'][combine_R_name]['three_prime_end_read2_adapter'] = df.loc[i, 'R_lib_R2_3prime_adaptor_seq']

        # extend adapter seq extract
        primer_F_seq = get_tag_seq(config['adapters'][combine_F_name]['five_prime_end_read2_adapter'], insert_tag, uplib=True)
        primer_R_seq = get_tag_seq(config['adapters'][combine_R_name]['five_prime_end_read2_adapter'], insert_tag, uplib=False)
        config['extend_seq'][sample] = [primer_F_seq, primer_R_seq]
    return config

def write_yaml(yaml_file, config):
    with open(yaml_file, 'w', encoding='utf-8') as f:
        yaml.dump(config, f)

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--data', type=str, default=None, required=True,
                        help='path for input file')

    parser.add_argument('-o', '--output', type=str, default=None, required=True,
                        help='path for output file')
    
    parser.add_argument('-k', '--kmer', type=int, default=15,
                        help='the read sequencing primer length for cutadapt, default 15bp')

    parser.add_argument('-t', '--thread', type=int, default=1,
                        help='the thread number for read mapping')

    args = parser.parse_args()

    data = args.data
    output = args.output
    # print(Path(output).resolve())
    Path(output).resolve().parent.mkdir(parents=True, exist_ok=True)

    df = read_table(data)
    config = generate_config(df, args)
    write_yaml(output, config)


if __name__ == '__main__':
    main()