# SubgraphMatching
## Introduction
Our ICDE'2024 paper "IVE: Accelerating Enumeration-based Subgraph Matching via Exploring Isolated Vertices".

If you have any further questions, please feel free to contact us.

Please cite our paper, if you use our source code.

* "IVE: Accelerating Enumeration-based Subgraph Matching via Exploring Isolated Vertices. ICDE 2024."


## Execution
Execute the following command to run the program.
```bash
make
./main -d ./data/data.txt -q ./data/query.txt -num 100000
```

## Input
Both the input query graph and data graph are vertex-labeled.
Each graph starts with 't N M' where N is the number of vertices and M is the number of edges. A vertex and an edge are formatted
as 'v VertexID LabelId Degree' and 'e VertexId VertexId' respectively. Note that we require that the vertex
id is started from 0 and the range is [0,N - 1] where V is the vertex set. The following
is an input sample. You can also find sample data sets and query sets under the test folder.

Example:

```bash
t 5 6
v 0 0 2
v 1 1 3
v 2 2 3
v 3 1 2
v 4 2 2
e 0 1
e 0 2
e 1 2
e 1 3
e 2 4
e 3 4
```

## Experiment Datasets
We use the [datasets](https://drive.google.com/file/d/1uHaVbNvZzUJipfw-tq-RH9lXeOTYD_hD/view?usp=sharing) in our experiments.


## References
[1] Zite Jiang, et al. Fast Subgraph Matching by Dynamic Graph Editing. TSC 2023.

[2] Shixuan Sun and Qiong Luo. In-Memory Subgraph Matching: an In-depth Study. SIGMOD 2020.
