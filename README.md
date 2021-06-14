# Bayesian Optimization for real-time, automatic design of face stimuli in human-centred research
[![MIT license](http://img.shields.io/badge/license-MIT-brightgreen.svg)](https://github.com/PedroFerreiradaCosta/FaceFitOpt/blob/master/LICENSE)


Official script for the paper ["Neuroadaptive electroencephalography: a proof-of-principle study in infants"](https://arxiv.org/abs/2106.06029).

![Paper Figure](./Figures/Figure1.png)

## Abstract
A core goal of functional neuroimaging is to study how the environment is processed in the brain. The mainstream paradigm involves concurrently measuring a broad spectrum of brain responses to a small set of environmental features preselected with reference to previous studies or a theoretical framework. As a complement, we invert this approach by allowing the investigator to record the modulation of a preselected brain response by a broad spectrum of environmental features. Our approach is optimal when theoretical frameworks or previous empirical data are impoverished. By using a prespecified closed-loop design, the approach addresses fundamental challenges of reproducibility and generalisability in brain research. These conditions are particularly acute when studying the developing brain, where our theories based on adult brain function may fundamentally misrepresent the topography of infant cognition and where there are substantial practical challenges to data acquisition. Our methodology employs machine learning to map modulation of a neural feature across a space of experimental stimuli. Our method collects, processes and analyses EEG brain data in real-time; and uses a neuro-adaptive Bayesian optimisation algorithm to adjust the stimulus presented depending on the prior samples of a given participant. Unsampled stimuli can be interpolated by fitting a Gaussian process regression along the dataset. We show that our method can automatically identify the face of the infant's mother through online recording of their Nc brain response to a face continuum. We can retrieve model statistics of individualised responses for each participant, opening the door for early identification of atypical development. This approach has substantial potential in infancy research and beyond for improving power and generalisability of mapping the individual cognitive topography of brain function.



## Citation
If you find this code useful for your research, please cite:

      @misc{dacosta2021neuroadaptive,
      title={Neuroadaptive electroencephalography: a proof-of-principle study in infants}, 
      author={Pedro F. da Costa and Rianne Haartsen and Elena Throm and Luke Mason and Anna Gui and Robert Leech and Emily J. H. Jones},
      year={2021},
      eprint={2106.06029},
      archivePrefix={arXiv},
      primaryClass={q-bio.NC}
      }
