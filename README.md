This is an anonymous repository submitted along with a paper submission for evaluation and reproducibility.

## Anchor-Based Architecture for Robust DNA Data Storage and Retrieval

Abstract: DNA-based data storage offers an unprecedented solution to the growing demand for long-term, high-density data storage. However, challenges in data recovery and error correction during the DNA readout process remain significant obstacles to its widespread adoption. In this paper, we propose a novel architecture for DNA data storage that introduces anchors—predefined sequences strategically embedded in the middle of the DNA strands encoding the data. These mid-strand anchors serve as reference points during data recovery, enabling more accurate alignment and reconstruction of data in the presence of errors. We demonstrate how this anchor-based architecture enhances the robustness of DNA data recovery, reduces error rates, and improves overall system efficiency. Through simulations and experimental analysis, we demonstrate that our approach enhances the reliability of data reconstruction, contributing to the advancement of DNA storage as a potential solution for next-generation archival systems.

This repository is forked from the [DNAStorageToolkit](https://github.com/DNAStorageToolkit/DNAStorageToolkit.git), with our changes to the encoding and reconstruction implemented on top of the pre-existing pipeline.

## Installation

To get started with our DNA storage pipeline toolkit, you need to install the following:
> * Python 3
> * G++
> * pyspoa
> * editdistance



After installing these, you can simply clone this repo and use the modules directly:

```
$ git clone https://github.com/DNA-Anchors/DNAnchors.git

$ cd DNAnchors
```
To see the pipeline in action, run the demo made available directly:
```
$ sh experiment.sh
```


## Usage

To demonstrate how to use this code, we have provided a demo script ([**experiment.sh**](./experiment.sh)) that takes the same image through the baseline DNA storage pipeline and our Anchored version of the pipeline.



## Citation

The baseline toolkit is borrowed from this work:
```
@inproceedings{sharma2024dnatoolkit,
  title={DNA Storage Toolkit: A Modular End-to-End DNA Data Storage Codec and Simulator},
  author={Sharma, Puru and Goh Yipeng, Gary and Gao, Bin and Ou, Longshen and Lin, Dehui and Sharma, Deepak and Jevdjic, Djordje},
  booktitle={2024 IEEE International Symposium on Performance Analysis of Systems and Software (ISPASS)},
  year={2024},
  abbr={ISPASS'24},
  organization={IEEE}
}
```
