# Maximum Balanced $(k, \epsilon)$-Bitruss Detection in Signed Bipartite Graph

Author : Qixiao ZHONG

Institution: Hong Kong University of Science and Technology

E-MAIL: qzhongaf@cse.ust.hk

# Usage
./blcbt datafile algorithm k epsilon

argv[0] = ./blcbt

argv[1] = datafile

	*   each line in format (uid vid sign)
    *   uid and vid are positive integers
    *   sign is either 1 or -1

argv[2] = algorithm

    *   countBL: baseline balanced butterfly counting algorithm.
    *   count: improved balanced butterfly counting algorithm (CountI).
    *   supportsGreedy: greedy heuristics by balanced supports ratio (GreedyS).
    *   followersGreedy: greedy heuristics by followers (GreedyF).
    *   exact: the exact algorithm.
	
argv[3] = $k$

    *   not required for countBL or count

	
argv[4] = $\epsilon$

    *   not required for countBL or count
	
Examples:

Run improved balanced butterfly counting algorithm

	*   ./blcbt datafile count

Run GreedyS algorithm with $k = 20$ and $\epsilon = 0.3$ 
	
	*   ./blcbt datafile supportsGreedy 20 0.3