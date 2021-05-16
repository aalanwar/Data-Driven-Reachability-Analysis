## Data Driven Reachability Analysis using Matrix Zonotopes
<br /> 
This repo cotains the code for our two papers:<br /> 
[1] Amr Alanwar, Anne Koch, Frank Allgöwer, Karl Johansson "Data Driven Reachability Analysis Using Matrix Zonotopes" 
3rd Annual Learning for Dynamics & Control Conference <br />
[2] Amr Alanwar, Anne Koch, Frank Allgöwer, Karl Johansson "Data Driven Reachability Analysis from Noisy" Submitted to IEEE Transactions on Automatic Control<br />


We consider the problem of reachability analysis from noisy data, given that the system model is unknown. 
Identifying a model is a preliminary step in the state-of-the-art reachability analysis approaches. 
However, systems are becoming more complex, and data is becoming more readily available. 
We propose a data-driven reachability analysis using matrix zonotope and using a new set representation named constrained matrix zonotope.<br />
The following figure summarize the idea behind our papers.
<br /> <br />
<p align="center">
<img
src="Figures/idea3.png"
raw=true
alt="Subject Pronouns"
width=500
/>
</p>
<br />
<br />
There are two levels of complexity for the proposed data driven reachability analysis<br />
A- Basic reachability analysis (files start with a-)<br />
B- Constrained matrix zonotope reachability analysis (files start with b-), 
These files compare three methods for reachability analysis namely, matrix zonotope, constrained matrix
zonotop using exact noise description and constrained matrix zonotope given side information.<br />
<br />

## Running 
1- Download [CORA](https://github.com/TUMcps/CORA) and [MPT](https://www.mpt3.org) toolboxs.<br />
2- Add CORA nad MPT folder and subfolders to the Matlab path.  <br />
3- Add the repo folder and subfolders to the Matlab path.  <br />
<br />
<br />
4- run a_linearDT.m for linear system using matrix zonotope.<br />
5- run a_nonlinearDT.m for nonlinear system.<br />
6- run a_polyDT.m for polynomial system using matrix zonotope.<br />
<br />
<br />
4- run b_linearDT_measnoise.m for linear system with measurement noise.<br />
5- run b_linearDT_sideInfo.m for linear system given side information.<br />
6- run b_polyDT_sideInfo.m for polynomial system given side information.<br />
<br />
<br />
You can save the workspace after any file starts with b and then run the plotting file which 
starts with p. For example, run<br />
b_linearDT_sideInfo.m<br />
save the workspace and load it later, run<br />
p_plot_linearDT_sideInfo.m<br />
<br />
<br />
<br />
Our papers Bibtex are as follow:<br />
```
@article{datadriven_reach1,
  title={Data-Driven Reachability Analysis Using Matrix Zonotopes},
  author={Alanwar, Amr and Koch, Anne and Allgöwer, Frank and Johansson, Karl Henrik},
  journal={arXiv preprint arXiv:2011.08472},
  year={2020}
}
```

```
@article{datadriven_reach2,
  title={Data-Driven Reachability Analysis from Noisy Data},
  author={Alanwar, Amr and Koch, Anne and Allgöwer, Frank and Johansson, Karl Henrik},
  journal={arXiv preprint},
  year={2021}
}
```