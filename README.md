# Temporal-k-Core-Query-Project

## ABSTRACT

Querying cohesive subgraphs on temporal graphs with time constraints is valuable since it reveals underlying communities that emerge during specific meaningful time intervals. In this paper, we study a general Temporal k-Core Query (TCQ) problem: find all distinct k-cores that exist within a given time interval in a temporal graph, which can be extended to many interesting applications. This problem is challenging due to the massive number of subintervals that possibly induce different k-cores. For that, we propose a novel Temporal Core Decomposition (TCD) algorithm that decrementally induces temporal k-cores from the previously induced ones and thus reduces **intra-core** redundant computation significantly. Then, we introduce an intuitive property named Tightest Time Interval (TTI) of temporal k-core, and design an optimization technique that leverages TTI to predict which subintervals will induce identical k-cores and safely prunes the subintervals in advance, thereby eliminating **inter-core** redundant computation. The complexity of optimized TCD (OTCD) algorithm no longer depends on the span of query time interval but only the scale of final results, which means OTCD algorithm is scalable. Moreover, we propose a compact in-memory data structure named Temporal Edge List (TEL) to implement OTCD algorithm efficiently in physical level. TEL can be updated instantly with new incoming temporal edges, and thus our approach can also deal with dynamic temporal graphs. We compare OTCD algorithm with the latest historical k-core query on several real-world temporal graphs, and observe that OTCD algorithm outperforms it by three orders of magnitude, while OTCD algorithm needs none precomputed index.

## File Description
1. TCD.cpp: Source code of TCD, OTCD
2. i-PHC.cpp: Source code of iPHC-TCQ

## Documentation
Experiments are conducted on Windows10 x64 with 64GB memory with Visual Studio 2019. At least 32GB memory is recommended. The code is implemented with pure C++ including C++17 features, so a compiler supporting latest C++ feature is recommended. The confirmed steps of compiling the source code are provided below:

1. Create a VS 2019 project and add the source code(i-PHC.cpp OR TCD.cpp) to the project.
2. Build the project with x64 configuration and generate the executable file.
3. Add the graph file and testcase file for testing to the directory of executable file.
4. Run the executable file with command.

The command for running the executable file is set as follows:```.\[exe file] [graph file] [testcase file]```, where graph file should only contain the edge set in form of "u v t", with **SPACE** as separator, and testcase file should only contain a group of queries line by line, in form of "ts te k" with **TAB** as separator. Please refer to the given example: graph-format.txt, query-format.txt.

After running the executable file, queries in the testcase file will be processed one by one through TCD and OTCD. For each processed query, multiple lines containing information of the query will be printed:
1. The name of graph
2. The processed query
3. The time consumed by TCD
4. The time consumed by OTCD
5. Total number of obtained temporal k-core in query

## Reproduction
**Hardware**: Windows10 x64 machine with 64 gigabytes memory and 20 gigabytes space.

To facilitate reproduction of experiments in paper, we create two packages which contains customized vs2019 project of iPHC-TCQ and TCD/OTCD repectively. The link is shared below:

**link**: https://pan.baidu.com/s/1cCG0jOofbi1jW1-2XmGxbw

**password**: dn3s

The default configuration is Debug with x64, under whose output file directory, we preserve all graphs and queries needed for reproducing corresponding figures or tables in paper. Recommended steps of reproduction are concluded below:

1. Build the project in Debug and x64.
2. Specify the graph file and query file by adding debugging arguments in project's attributes.
3. Press **"local windows debugger"** button to start running.

For example, to reproduce the diagram of graph CollegeMsg in Figure 7 of paper, set the debugging arguments as: 

``.\Figure7\graph\CollegeMsg.txt .\Figure7\query\CollegeMsg-q.txt``.

Macros are also adopted for saving time and memory in reproduction. When reproducing Table 4, enable the preprocessor's option ``CALC_PRUNE`` than rebuild the project, otherwise disable it.

Output information for each query without ``CALC_PRUNE`` are:

1. Graph Name: The graph name
2. Platform: (Ignore this)
3. Query: The processed query [ts, te, k]
4. TCD Time Clapse: Time consumption of TCD for resolving the query
5. TCD Accessed Cell: Number of Cell accessed by TCD
6. OTCD Time Clapse: Time consumption of OTCD for resolving the query
7. OTCD Accessed Cell: Number of Cell accessed by OTCD

If ``CALC_PRUNE`` is enabled, a few extra lines will be included:

1. PoR: Times of PoR triggered
2. PoU: Times of PoU triggered
3. PoL: Times of PoL triggered
4. CPoR: Number of cells pruned by PoR
5. CPoU: Number of cells pruned by PoU
6. CPoL: Number of cells pruned by PoL

## Contact
If you have any questions, contact us by sending an email to clock@whu.edu.cn / thomasyang@whu.edu.cn
