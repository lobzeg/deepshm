# DeepSHM
This repository implements the model from "Deep learning model of somatic hypermutation reveals importance of sequence context beyond targeting of AID and Polη hotspots" by Tang, Krantsevich & MacCarthy (https://www.biorxiv.org/content/10.1101/2021.08.03.453264v1.full.pdf). Using DeepSHM you can simply use one of the provided models to make SHM-related predictions for your data, or you can train your own CNN model for any task that uses DNA/RNA one-hot encoded data as input.


## Instalation
### If you want to play with the model 
DeepSHM is on pypi, so it can be installed using pip:

    pip install deepshm

### If you are not familiar with python
You can also just download the deepshm_predict.zip archive and 0nce you extract it you can call the file deepshm_predict.py from the command line. Note that you still will need python installed as wells as several packages: tensorflow, numpy, scipy. deepshm_predict.py can be called in the following way:

    python deepshm_predict.py input.fasta -o output.csv -k 15

where the input fasta file is the required input and the path to the output file and the length of k-mer (5, 9, 15 and 21) are optional. Default value for the output file is output.csv and default k is 15.
