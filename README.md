# PathWalker

In this work, we present a robust algorithm, PathWalker, that improves the power of RWR method in reconstructing pathways by inspiring ideas from the PathLinker algorithm, resulting in a pathway with an accuracy higher than both algorithms.

The pathway constructed by the RWR algorithm is concentrated over a neighborhood of the set of receptors and up until the first few hundred edges, does not find paths connecting receptors to the set of transcription factors. The reason for this phenomenon is that in the RWR algorithm, the choice of picking the next vertex is independent of the set of transcription factors; i.e, if the random walker is at a vertex u, then he/she chooses the next vertex to visit randomly from the set of neighbors of u. Instead, in the PathWalker algorithm, we give a higher priority to the neighbors of u that are closer to the set of transcription factors. Next, we describe the details of the algorithm.

## Step 1: Computing the distances
In the first step of the algorithm, we want to compute the probability of the highest-probability path from each node u to any transcription factor in the set of transcription factors T. To do so, we construct an augmented interactome in which, 
1. The edge weights are negative log_2 transformed; i.e, for each edge (u,v) with weight Wuv, we replace the edge weight with -log2 Wuv , 
2. The direction of all edges are reversed, and  
3. A fake source vertex s with zero-weight edges connecting s to the vertices in T is added.
   
We execute Dijkstra's shortest path algorithm from the fake source s. For each vertex v, suppose lv denotes the length of the shortest path returned by Dijkstra's algorithm. Define Sv is the probability of the highest-probability path from each node v to the set of transcription factors in T.
## Step 2: Adjusting the edge weights
In this step, we redefine the edge weights. For any edge directed from a node u to a node v, Define Suv = Su Puv to be the probability of the highest-probability path connecting u to any transcription factor that passes through the node v. We set the weight of the edge (u,v) to Wuv = Suvr , where r is a parameter determining the importance of the higher-probability paths; i.e, for each node u, considering that the random walker is at node u, the higher values for r makes the random walker choose the edges (u,v) with higher values of Suv with a higher probability. 
## Step 3: Running the RWR algorithm with the adjusted edge weights 
Finally, we simply run the random walk with restarts algorithm on the network with edge weights defined as above to obtain the stationary probability Pv for each node v of the network, followed by assigning the priorities and sorting the edges based on their priorities.

