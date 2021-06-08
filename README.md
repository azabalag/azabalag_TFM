# TFM_UOC

Ésto es una guía de los códigos empleados en el Trabajo de Fin de Máster titulado "Desarrollo y mejora de la funcionalidad de la herramienta _XICRA_: comparativa de análisis de _tRNAs_"

## XICRA

### Instalación 

La instalación de _XICRA_ (https://github.com/HCGB-IGTP/XICRA) ha de hacerse en modo desarrollo:

```bash
## clone repo
git clone https://github.com/HCGB-IGTP/XICRA.git

## move to folder XICRA_pip
cd XICRA/XICRA_pip

## create conda environemt
conda env create -n XICRA -f ./devel/conda/environment.yml

## activate
conda activate XICRA

## install latest python code
pip install -r ./devel/pypi/requirements.txt
pip install -e .

## install missing software
sh ./XICRA/config/software/installer.sh
```

Para chequear si se ha realizado adecuadamente la instalción:

```bash
XICRA config
```

### Análisis de _miRNA_

Se adjunta el código de ejemplo de uso del análisis de _miRNA_:

```bash
## run XICRA example
ln -s ~/BMC_bioinformatics_paper/simulation/example/reads/

## prepare reads
XICRA prep --input reads/ --output_folder test_XICRA

## join reads
XICRA join --input test_XICRA --noTrim

## create miRNA analysis
XICRA miRNA --input test_XICRA --software miraligner sRNAbench

## explore results
ls test_XICRA/report/
```

## Obtención de datos:

Se detalla el código para obtener las muestras con _SRAtoolkit_ (https://github.com/ncbi/sra-tools). Se muestra como ejemplo la obtención de la muestra SRR12344552:

```bash
fasterq-dump --split-files SRR12344552
```

### Preparación de los datos

Es necesario adecuar las muestras descargadas para realizar el análisis de estos. Estos pasos se realizaron con los módulos de XICRA:

```bash
## prepare reads
XICRA prep --input reads/ --output_folder test_XICRA

## trimming
XICRA trimm --input --adapters_a TGGAATTCTCGGGTGCCAAGG --adapters_A GATCGTCGGACTGTAGAACTCTGAAC

## join reads
XICRA join --input test_XICRA
```

## Instalación de los softwares
 
### _MINTmap_

La instalacion de _MINTmap_ (https://github.com/TJU-CMC-Org/MINTmap) se realizó con Conda:

```bash
conda install -c bioconda mintmap
```

### _miRge3.0_

La instalación de _miRGe3.0_ (https://github.com/mhalushka/miRge3.0) se tuvo que realizar en un entrono de Conda propio:

```bash
## create conda environemt
conda env create -n mirge

## activate
conda activate mirge
```

```
## install python3.8 and R
sudo apt update
sudo apt install software-properties-common
sudo add-apt-repository ppa:deadsnakes/ppa
sudo apt install python3.8
sudo apt install python3-setuptools
sudo apt install python3-pip
sudo apt install r-base
```

```bash
conda install -c bioconda mirge
```

Así como su dependencias:

```bash
## bowtie
wget -O bowtie-1.3.0-linux-x86_64.zip https://sourceforge.net/projects/bowtie-bio/files/bowtie/1.3.0/bowtie-1.3.0-linux-x86_64.zip/download
unzip bowtie-1.3.0-linux-x86_64.zip
cd bowtie-1.3.0-linux-x86_64
pwd 
  /home/arun/software/bowtie-1.3.0-linux-x86_64
export PATH=$PATH:"/home/arun/software/bowtie-1.3.0-linux-x86_64"

## samtools
sudo apt install samtools

## RNAfold
wget “https://www.tbi.univie.ac.at/RNA/download/sourcecode/2_4_x/ViennaRNA-2.4.16.tar.gz”
cd ViennaRNA-2.4.16
sudo ./configure 
sudo make 
sudo make install

## GUI requirements
sudo ln -s /home/arun/.local/bin/miRge3.0 /usr/local/bin/miRge3.0
sudo ln -s /home/arun/.local/bin/cutadapt /usr/local/bin/cutadapt
sudo ln -s /home/arun/software/bowtie-1.3.0-linux-x86_64/bowtie /usr/local/bin/bowtie
sudo ln -s /home/arun/software/bowtie-1.3.0-linux-x86_64/bowtie-build /usr/local/bin/bowtie-build
sudo ln -s /home/arun/software/bowtie-1.3.0-linux-x86_64/bowtie-inspect /usr/local/bin/bowtie-inspect
```

### _tDRmapper_

La instalción de _tDRmapper_ (https://github.com/sararselitsky/tDRmapper) se realizo con Conda:

```bash
conda install -c bioconda tdrmapper
```

## Analisis de _tRNA_

Se muestra como ejemplo el análisis de _tRNA_ para la muestra sRR12344552

### _MINTmap_

```bash
MINTmap.pl -f .SRR12344552_trim_joined.fastq -l ./LookupTable.tRFs.MINTmap_v1.txt -s ./tRNAspace.Spliced.Sequences.MINTmap_v1.fa -o ./OtherAnnotations.MINTmap_v1.txt 
```

### _miRGe3.0_

```bash
./miRge3.0 -s /./SRR12344552_trim_joined.fastq -lib ./miRge3_Lib/ -on human -db mirgenedb -pbwt -trf -gff -ie
```
### _tDRmapper_

```bash
perl /./TdrMappingScripts.pl ./hg19_mature_and_pre.fa ./SRR12344552_trim_joined.fastq 
```

## Comparación de la detección de _tRFs_

A continuación se detalla el usos de los scripts creados para realizar la comparación de la detección de _tRFs_. Como ejemplo se muestra la comparación de la muestra SRR12344552:

### table.py

Resume los resultados obtenidos por cada software en una tabla, y muestra las secuencias _tRFs_ detectadas y los contajes crudos detectados de cada uno de ellos.

```bash
## MINTmap Ambiguous
python table.py MINTmapAmbiguous SRR12344552 /home/zabala/Escritorio/Test_XICRA/MINTmap_SRR12344552/output-MINTmap_v1-ambiguous-tRFs.expression.txt

## MINTmap Exclusive
python table.py MINTmapExclusive SRR12344552 /home/zabala/Escritorio/Test_XICRA/MINTmap_SRR12344552/output-MINTmap_v1-exclusive-tRFs.expression.txt

## miRGe3.0
python table.py miRge3.0 SRR12344552 /home/zabala/Escritorio/Test_XICRA/miRge/tRFs.samples.tmp/SRR12344552_trim_joined.tRFs.report.tsv 

## tDRmapper
python table.py tDRmapper SRR12344552 /home/zabala/Escritorio/Test_XICRA/tDRmapper/SRR12344552_trim_joined.fastq.hq_cs.mapped
```

### summary.py

Resume los _tRFs_ totales y únicos del software analizado.

```bash
python summary.py SRR12344552
```

### upsetR.R

Compara los _tRFs_ detectados en común por los softwares, y los muestra en un gráfico:

```r
library(UpSetR)

upsetR(MINTmap_Ambiguos, MINTmap_Exclusive, miRGe3.0, tDRmapper)
```

### comparation.py

Detecta los tRFs detectados en común por los softwares. Esto es, cuantifica el solapamiento de los softwares en la detección de tRFs.

```bash
conda comparation.py SRR12344552
```
