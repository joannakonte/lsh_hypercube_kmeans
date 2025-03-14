# Approximate k-NN & Clustering on MNIST  

## Overview  
This project explores **approximate k-nearest neighbors (k-NN)** and **Approximate Range Search** using:  
- **Locality-Sensitive Hashing (LSH)**  
- **Randomized Projection on Hypercubes**  
- **K-Means++ Clustering**  

The **MNIST dataset** of handwritten digits is used for evaluation, with Euclidean distance as the primary metric. The project is implemented in **C++**.

It was done in 2023 as an assignement in the course "Software Development for Algorithmic Problems" by Panagiotis Chatzimichos and Ioanna Kontemeniotou. 

## Repository Structure  
```
/build              # Compiled binaries  
/include            # Header files  
/modules            # Core implementation (LSH, clustering, hypercube)  
/src                # Main program files  
Makefile            # Build system  
cluster.conf        # Configuration file  
report.pdf          # Detailed project report  
```  

## Compilation & Execution  
```sh
make                # Build the project  
make run_lsh        # Run LSH  
make run_cube       # Run Hypercube  
make run_cluster    # Run Clustering  
make clean          # Remove compiled files  
```  

## Parameters  
Modify values in the Makefile or use command-line arguments:  
- **LSH:** `-d <input_file> -q <query_file> -k <hash_functions> -L <hash_tables>`  
- **Hypercube:** `-d <input_file> -q <query_file> -k <projection_dim> -M <max_candidates>`  
- **Clustering:** `-i <input_file> -c <config_file> -o <output_file> -m <method>`  

## More Information  
For a detailed explanation of the implementation, experiments, and results, please refer to **[report.pdf](report.pdf)**.  

## Contributors 

• [Panagiotis Chatzimichos](https://github.com/pchatz000)\
• [Ioanna Kontemeniotou](https://github.com/joannakonte)