## CatMOD: a CatBoost-ensemble decision tree framework designed for  m<sup>6</sup>A modification recognition from ONT DRS

![Static Badge](https://img.shields.io/badge/Linux-blue?logo=Linux&logoColor=white)![](https://camo.githubusercontent.com/cd1d74df39562c27a5d2700b62a3c29bf1fc8469f1485c36fba8a06336604a3f/68747470733a2f2f696d672e736869656c64732e696f2f62616467652f6c616e67756167652d707974686f6e2d626c75652e737667)

## CatMOD overview

<b>CatMOD</b>, a CatBoost-ensemble decision tree framework designed for single-base resolution detection of m<sup>6</sup>A modifications , trained on a dataset supported by DRS and m<sup>6</sup>A -seq profiles. The model incorporates four feature dimensions: <b>biological genomic sequences</b>, <b>base-calling qualities</b>, <b>systematic alignment errors</b>, and <b>current signals</b>.
![Static Badge](https://github.com/cma2015/CatMOD/blob/main/img/catmod_workflow.png)

## Run CatMOD locally
### Requirements
CatMOD Project is a python3 package. To use CatMOD, python version 3.9 or higher is required.

- python >= 3.9
- catboost
- h5py
- numpy
- pysam
- rich
- scipy

###  Installation
```
git clone https://github.com/CatMOD/CatMOD.git
cd CatMOD
conda create -n catmod -y python=3.9
python setup.py install

##########
git clone https://github.com/CatMOD/CatMOD.git
cd CatMOD
conda env create -f catmod.yml

##########
git clone https://github.com/CatMOD/CatMOD.git
conda create -n catmod -y python=3.9 catboost h5py numpy pysam rich scipy
conda activate catmod
```
## How to use CatMOD
### PreProcess
```
### basecalling
guppy_basecaller --input_path $fast5_folder --recursive --fast5_out --save_path $guppy_folder --flowcell $FLOWCELL --kit $KIT --num_callers $THREADS

### multi-fast5 to single-fast5
multi_to_single_fast5 --input_path $guppy_folder --save_path $single_folder --threads $THREADS --recursive

### tombo resquiggling
tombo resquiggle --rna --processes $threads --overwrite --fit-global-scale --include-event-stdev $single_folder $REFERENCE
```

### Extracting Features
```
catmod extract_features --bed $sample_bed --ref $REFERENCE --align $ont_bam --current $ont_current --threads $THREADS --output $datasets_folder
```

### Predicting
```
catmod predict --bed $sample_bed --datasets $datasets_folder --model /path/to/CatMOD/models/wheat_pretrained.cbc.cbm --threads $THREADS --output $datasets_folder
```

## How to access help
- Comments/suggestions/bugs/issues are welcome reported [here](https://github.com/cma2015/CatMOD/issues) or contact: Minggui Song smg@nwafu.edu.cn or Chuang Ma chuangma2006@gmail.com

## Change log
- 2025.04 Release CatMOD v1.0
- 2022.06 we launched CatMOD project
