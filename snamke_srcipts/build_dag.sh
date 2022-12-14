#########################################################################
# File Name: build_dag.sh
# Author: luok
# mail: luokai@westlakegenetech.com
# Created Time: 2022年11月15日 星期二 13时37分36秒
#########################################################################
#!/bin/bash

# snakemake -p -s ../offtarget_pipline.sh -j 24 -c 64

snakemake  --dag -s offtarget_pipline.sh   |dot -Tpdf > pipeline.pdf
