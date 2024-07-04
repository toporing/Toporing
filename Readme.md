# Toporing

Repository for the source code of the engine presented in the paper Worst-Case-Optimal Joins on Graphs with Topological Relations 

## Instructions

To run our code, **we have to install an extended version of the library SDSL**. Go to [this repository](https://github.com/toporing/sdsl-lite) and follow the instructions.

After the extended version of SDSL is installed, we have to clone this repository and follow these steps:

1. Create our `build` folder and compile the code:
```Bash
mkdir build
cd build
cmake ..
make
```

Check that there is no errors.

2. Download the version of Wikidata used in the testing:

- [Wikidata Filtered (957,450,487 triples)](http://doi.org/10.6084/m9.figshare.25683684).

Now put the .dat file inside a folder.

3. Building the index. After compiling the code we should have an executable called `build-index` in `build`. Now run:

```Bash
./build-index <absolute-path-to-the-.dat-file> <CONTAINS_ID> <CONTAINED_IN_ID> <TOUCHES_ID>
```


4. Querying the index. In `build` folder, you should find another executable file called `query-index`. To solve the queries you should run:

```Bash
./query-index <absoulute-path-to-the-index-file> <absolute-path-to-the-query-file> <LIMIT_RESULTS> [BASELINE]
```

By setting <LIMIT_RESULTS> to 0, queries are executed until completed or timed_out. [BASELINE] can be set to 1 to run our baseline. This option is turned off by default 

Note that the second argument is the path to a file that contains all the queries. The queries of our benchmark are in `Queries`:

After running that command, you should see the number of the query, the number of results, and the elapsed time of each one of the queries with the following format:
```Bash
<query number>;<number of results>;<elapsed time>;<user time>
```
---
## Baselines

For running the baselines described in the paper, see [here](https://compact-leapfrog.imfd.cl/) for instructions.

