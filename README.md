# MLP simlator (+NeuroSim) V3.0



The MLP+NeuroSim framework was developed by [Prof. Shimeng Yu's group](http://faculty.engineering.asu.edu/shimengyu/) (Georgia Institute of Technology). The model is made publicly available on a non-commercial basis. Copyright of the model is maintained by the developers, and the model is distributed under the terms of the [Creative Commons Attribution-NonCommercial 4.0 International Public License](http://creativecommons.org/licenses/by-nc/4.0/legalcode)

This is the released version 3.0 (Feb, 2019) for the tool. In this version, the algorithm weight range is changed from (0,1) in V2.0 to (-1,1). Parrallel read-out architecture is added for digital eNVM based synaptic arrays. More training algorithms are available.

Developers: Pai-Yu Chen, Xiaochen Peng and Yandong Luo
PI: Prof. Shimeng Yu, Georgia Institute of Technology

If you have logistic questions or comments on the model, please contact Prof. Shimeng Yu (shimeng.yu@ece.gatech.edu), and if you have technical questions or comments, please contact Xiaochen Peng (xpeng76@gatech.edu) or Yandong Luo (yluo310@gatech.edu)

This research is supported by NSF CAREER award, NSF/SRC E2CDA program, and ASCENT, one of the SRC/DARPA JUMP centers. 

If you use the tool or adapt the tool in your work or publication, you are required to cite the following reference:

P.-Y. Chen, X. Peng, S. Yu, ※NeuroSim+: An integrated device-to-algorithm framework for benchmarking synaptic devices and array architectures,*§ IEEE International Electron Devices Meeting (IEDM)*, 2017, San Francisco, USA.

## File lists
1. MATLAB fitting script: `nonlinear_fit.m`
2. Nonlinearity-to-A table: `Documents/Nonlinearity-NormA.htm`
3. MNIST data: `MNIST_data.zip`
4. Manual: `Documents/Manual.pdf`
5. MLP Simulator (+NeuroSim): the rest of the files

## Installation steps (Linux)
1. Get the tool from GitHub
```
git clone https://github.com/neurosim/MLP_NeuroSim_2.0.git
```

2. Extract `MNIST_data.zip` to it’s current directory
```
unzip MNIST_data.zip
```

3. Compile the codes
```
make
```

For the usage of this tool, please refer to the manual.

## References related to this tool 
1. P.-Y. Chen, X. Peng, S. Yu, ※NeuroSim+: An Integrated Device-to-Algorithm Framework for Benchmarking Synaptic Devices and Array Architectures,*§ IEEE International Electron Devices Meeting (IEDM), 2017, San Francisco, USA.
2. P.-Y. Chen, X. Peng, S. Yu, ※System-level benchmark of synaptic device characteristics for neuro-inspired computing,§ IEEE SOI-3D-Subthreshold Microelectronics Technology Unified Conference (S3S) 2017, San Francisco, USA.
3. P.-Y. Chen, S. Yu, ※Partition SRAM and RRAM based synaptic arrays for neuro-inspired computing,*§ IEEE International Symposium on Circuits and Systems (ISCAS)*, 2016, Montreal, Canada.
4. P.-Y. Chen, B. Lin, I.-T. Wang, T.-H. Hou, J. Ye, S. Vrudhula, J.-S. Seo, Y. Cao, and S. Yu, ※Mitigating effects of non-ideal synaptic device characteristics for on-chip learning,*§ IEEE/ACM International Conference on Computer-Aided Design (ICCAD)*, 2015, Austin, TX, USA.
5. P.-Y. Chen, D. Kadetotad, Z. Xu, A. Mohanty, B. Lin, J. Ye, S. Vrudhula, J.-S. Seo, Y. Cao, S. Yu, ※Technology-design co-optimization of resistive cross-point array for accelerating learning algorithms on chip,*§ IEEE Design, Automation & Test in Europe (DATE)*, 2015, Grenoble, France.stic questions or comments on the model, please contact [Prof. Shimeng Yu](mailto:shimengy@asu.edu), and if you have technical questions or comments, please contact [Pai-Yu Chen](mailto:pchen72@asu.edu) and [Xiaochen Peng](mailto:xpeng15@asu.edu).
