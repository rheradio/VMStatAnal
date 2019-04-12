=============================================================================
SUPPORTING THE STATISTICAL ANALYSIS OF VARIABILITY MODELS
=============================================================================
R. Heradio, David Fernandez-Amoros, Christoph Mayr-Dorn, and Alexander Egyed.
In 41st ACM/IEEE International Conference on Software Engineering (ICSE), 
Montréal, QC, Canada, 2019.
=============================================================================

=============================================================================
- Structure of the artifact ------------------------------------------------- 
=============================================================================

This artifact is structured in the following directories:

1) installer: precompiled binaries to facilitate the instalation of our software 
   in Linux and MacOS (see INSTALL.txt).
                   
2) source: includes the source code of our software. It is organized in two 
   subfolders:
   2.1) cudd-release: as noted in the paper, our algorithms have are implemented 
        as an extension of the library CUDD 3.0 for BDDs. This directory 
        includes such CUDD release.
   2.2) paper-src: source code for the Feature Inclusion Probability (FIP) and
       Product Distribution (PD) algorithms described in the paper.

3) models: contains the variability models in their original 
   notation (i.e., the paper experimental data):
   3.1) KConfig for:
     3.1.1) axTLS 1.5.3 (http://axtls.sourceforge.net/)
     3.1.2) Fiasco 2014092821 (https://os.inf.tu-dresden.de/fiasco/)
     3.1.3) uClibc 20150420 (https://www.uclibc.org/)
     3.1.4) Busybox 1.23.2 (https://busybox.net/)
     3.1.5) EmbToolkit 1.7.0 (https://www.embtoolkit.org/)
   3.2) The SPLOT feature modelling language for:
     3.2.1) Dell laptops (A. Nöhrer and A. Egyed, C2O configurator: a tool for 
            guided decision-making, Automated Software Engineering, vol. 20, 
            no. 2, pp. 265–296, 2013).
     3.2.3) Automotive 02 (S. Krieter, T. Thüm, S. Schulze, R. Schröter, and 
            G. Saake, Propagating configuration decisions with modal implication 
            graphs, ICSE'18, New York, NY, USA, 2018, pp. 898–909).
     
4) results: contains a subdirectory for each variability model tested. In each 
   of these subdirectories, you may find the following files:
   4.1) name.dddmp file: The BDD corresponding to the variabiblity model in 
        dddmp ASCII format. There is a header that contains relevant information 
        such as the number of nodes in the BDD and the names and ordering of 
        the variables.
   4.2) name.var: The order of the variables in the BDD.
   4.3) name.exp: The translation of the model to propositional logic. The dddmp
        file was generated from the .var and the .exp file by successively 
        "applying" the constraints with a BDD engine.
   4.4) name.probability: The results of computing the probability of each 
        feature appearing in a software product generated from the variability 
        model, i.e., the results of applying the FIP algorithm described in the 
        paper.
   4.5) name.histogram: The histogram of the variability model, i.e. how many 
        products of each size there are, i.e., the results of applying the
        PD algorithm described in the paper.
   4.6) histogram.log: A log of the execution of the histogram tool. It is only 
        relevant to know the running time.
   4.7) probability.log: A log of the execution of the probability tool (for 
        running times purposes).

=============================================================================
- How to use the artifact? -------------------------------------------------- 
=============================================================================

The usage is extremely simple: Just type the name of the tool (probability or
histogram), the name of the .dddmp file that stores the BDD encoding of a 
variability model (without the extension), and use redirection to get a results 
file. For instance:

probability axtls >output.probability

You may also use a multithreaded version by specifying the number of threads 
as follows:

probability -t 20 >output.20.probability

This would use 20 threads to process the BDD.

The results in the paper were obtained using the single-threaded options except 
for the embtoolkit variability model, for which 50 threads were employed. 
This is because whereas the topology of our embtoolkit BDD supports an efficient
parallelization, for the other BDDs a single thread execution has better 
performance than a parallel one.  
