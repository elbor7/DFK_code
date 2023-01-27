<h1>Discrete Polyhedral Cell Complex (PCC) Fracture Kinetic code (DFK code)</h1>

Manual version: 0.1.0 <br>
First manual release date: 21/09/2022 <br>
Current manual release date: 27 January 2023 <br>

<p> DFK code is a <i> fracture simulation tool </i> designed for polycrystalline composites design purposes. Its key feature is the usage of polyhedral cell complexes (PCC), which provide a discrete space for designing realistic material defect structures of different dimensions and types. Several functions allow the simulation physically realistic intragranular fracture process, including multiple cracking and macrocrack growth. Specifically, it was tested on the problems of complex-structured rGO/ceramic nanocomposites fracture.
</p>

<p> This is a C++ based software project consisting of several modules working with a pre-created Polyhedral Cell Complex (PCC) as the set of its incidence and adjacency matrices represented in a sparse matrix form. The code contains several modules (Processing, Subcomplex, Multiphysics, Kinetic, Characterisation and Writer) with additional libraries containing special classes and support functions. Each specific simulation problem (studies of the effect of elastic stresses, inclusions fraction, etc.) can be included in the main.cpp as a separate simulation task. The main.cpp merge and launch all the modules together, reading configurations and data files, "include" simulation tasks.

The code works both with 3-complexes (3D) and 2-complexes (2D) BUT - to make results consistent with 2D/3D EBSD scans - it assumes that the "grains" are 3-cells in the 3D, 2-cells in the 2D case, and so on for grain boundaries and other element types. So it actually replaces definitions of k-cells with ( k + (dim - 3) )-cells, where dim = {2, 3} is the dimension of the problem. In such a way, 2-cells are ALWAYS associated with grain boundaries and are edges for the 2D case.
</p>
  
<h2>General specifications</h2>
  
<p>The code is written and tested in C++ 17 with the parallelised verson used some features of C++ 20. It explicitly uses the Eigen and Spectra libraries, which must be <a href="https://spectralib.org/download.html"> downloaded</a> and copied to the directory containing all the STL C++ libraries on the local PC. The code works well and has been tested using CMake 3.23 (cmake.org), g++ compiler (gcc.gnu.org) and CLion IDE (jetbrains.com/clion).

The computational costs of different calculation types, functions and tasks are hugely different: for instance, the component analysis or the Metropolis algorithm-based simulations are very time-consuming procedures, while the random generation of special chains is fast.<\p>
  
<h2>Modules</h2>
All the modules except the main.cpp consist of the <i>interface<\i> and the core <i>library<\i> parts. The interfaces contain pre- and post-processing of data for this particular module, adapting input from the main.cpp, and managing the function implementations according to the calculation types specified in the configuration file. 

<ol>
  <li>DCC_Main: include libraries, global variables, reading from files, launching the other modules and simulation tasks; </li>
  <li>DCC_Processing: assign chains of special k-cells. Essentially, the ultimate goal of the module is to create a list of special k-cells (or k-sequence) in the order of appearance of new special elements (<i>inclusions</i>) during the consideration process; </li>
  <li>DCC_Subcomplex: including Plane cut (a,b,c,D), (reduced (k-1)-complex)</li>
  <li>DCC_Multiphysics: set energies for all the complex elements and especially its faces or 2-cells </li>
  <li>DCC_Kinetic: contain several functions implementing the process of multiple cracking. Essentially, the ultimate goal of the module is to create another list of special k-cells (or k-sequence) in the order of appearance of new <i>fractured</i> elements during the fracture kinetic process;   </li>
  <li>DCC_Characterisation: performs characterisation of special structures; </li>
  <li>DCC_Writer: formatted output simulations data in *.txt files; </li>
  <li>Objects.h: contain a library of C++ classes like the grain3D, subcomplex, macrocrack, etc. </li>
  <li>Support_Functions.h: contains all the other necessary support functions for the code </li> 
</ol>

<h3> Input files </h3>
 <p>
All the input files must be in a single folder specified as the ‘input’ directory in the ‘config.txt’ file.

<i>config.txt</i> file contains
<ol>
<li> ‘input’ directory address </li>
<li>  ‘output’ directory address </li>
<li>   the number of calculation points for different values of <i>p</i></li> and other parameters necessary for the code</li>
 </p>


<h2> More sources on the mathematics related to the Polyhedral cell complexes</h2>
An excellent simple introduction to the DCC with their various applications is given in the <a href="https://link.springer.com/book/10.1007/978-1-84996-290-2" target="_blank"> book </a> of Leo Grady and Jonathan Polimeni _“Discrete Calculus. Applied Analysis on Graphs for Computational Science. (2010)_

<h2> Where to take a complex? </h2>
The discrete cell complex is a pretty well-known object that originated from the field of algebraic topology, so it can be obtained in many various ways Below is just a concise review of a couple of flexible tools developed in the Mechanics and Physics of Solids research group in the University of Manchester providing DCCs based on Voronoi and a few others tessellations of space by convex polygons. 

<h3> Tessellations of space provided by Neper software </h3>
The Voronoi tesselation provided by Neper supposed to be a <i>dual</i> complex and so all the other tessellations provided by the Neper output with the <a href="https://neper.info/doc/neper_t.html#morphology-options" target="_blank"> morphology </a> option <i> -morpho <morphology> </i> like <i> cube, square, tocta, lamellar, etc. </i> different from <i>voronoi</i>.

Please, see more <a href="https://neper.info/doc/neper_t.html#examples" target="_blank"> examples </a> on the Neper webpage.

<h2> Applications of the DСD tool </h2>
<ol>
<li> E.N. Borodin, A.P. Jivkov, A.G. Sheinerman, M.Yu. Gutkin, 2021. Optimisation of rGO-enriched nanoceramics by combinatorial analysis. Materials & Design 212, 110191. [doi: 10.1016/j.matdes.2021.110191.](https://doi.org/10.1016/j.matdes.2021.110191) </li>
</ol>

<h2> Acknowledgements </h2>
This code was created by Dr Elijah Borodin in 2022 for the project RSF 18-19-00255 (https://rscf.ru/project/18-19-00255/)

<h2> License </h2>
The code is distributed under the GNU General Public License v3.0. See LICENSE.txt for more information.

<h2> Contacts </h2>
Dr Elijah Borodin (Research Fellow in Materials Physics at the University of Manchester; Mechanics and Physics of Solids research group)
<a href=“ Elijah.Borodin@icloud.com” Send e-mail> to Elijah Borodin (any queries regarding the code) 
