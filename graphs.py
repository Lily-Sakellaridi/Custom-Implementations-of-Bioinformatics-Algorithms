#!/usr/bin/env python3

"""
Author: Lily Sakellaridi
Given a directed, connected graph G, the script:
1) determines if G is eulerian
2) determines if G contains an eulerian path
3) finds an eulerian cycle in G, assuming one exists
4) finds an eulerian path in G, assuming one exists.

Given a spectrum, the script constructs the graph that
corresponds to the spectrum. Then it determines if that graph is
eulerian or has an eulerian path, finds the eulerian path,
and prints the sequence it corresponds to.

The algorithm and referenced figures can be found in Chapter 8 of 
"An Introduction to Bioinformatics algorithms" (Jones & Pevzner).
"""
# Import statements
from sys import argv, exit
from itertools import permutations
from random import choice

# Implement your functions here

def find_all_nodes(graph):
    """
    Return all nodes in a graph.
    
    graph: directed dictionary
    keys: nodes; strings or integers
    values: neighbor nodes; lists of strings or lists of integers
    """
    nodes = []
    
    for node, neighbors in graph.items():
        if node not in nodes:
            nodes.append(node)
        for neighbor in neighbors:
            if neighbor not in nodes:
                nodes.append(neighbor)
                
    return nodes
        
def is_eulerian(graph):
    """
    Return True if the provided graph is Eulerian, False otherwise.
    
    graph: directed graph; dictionary where:
    keys: nodes; strings or integers
    values: lists of connected nodes (strings or integers)
    
    Explanation: for a graph to be Eulerian, two requirements need
    to be fulfilled:
    1) The graph is connected
    2) Every node of the graph is balanced.
    
    The function assumes that the provided graph is connected,
    and checks only the second requirement.
    
    To show that a node is balanced, outgoing and incoming edges are
    counted. If they are equal in number for all nodes, the graph
    is Eulerian.
    """
    
        
    nodes = find_all_nodes(graph)
    eulerian = True
    
    for node in nodes:
        outdegree = len(graph.get(node, "")) # evaluates to 0 if node has only incoming edges
        indegree = len([value for value in graph.values() if node in value])
        if outdegree != indegree:
            eulerian = False
            break
            
    return eulerian
    
def has_eulerian_path(graph):
    """
    Return True if graph has an Eulerian path, False otherwise.
    
    graph: directed graph; dictionary
    keys: nodes; strings or integers
    values: neighbor nodes; lists of strings or lists of integers
    
    A graph has an Eulerian path is it contains at most two (2)
    semi-balanced nodes.
    'Semi-balanced' means that abs(indegree - outdegree) == 1
    
    The function iterates over the nodes, checking
    the indegree and outdegree of each.
    
    If they are equal, the node is balanced and therefore ignored.
    If not, the node is checked to see if it's semi-balanced. A 
    counter for semi-balanced nodes is kept.
    
    The function returns False in two cases:
    1) a node is found that is neither balanced nor semi-balanced
    2) the number of semi-balanced nodes becomes larger than 2.
    """
    nodes = find_all_nodes(graph)
    eulerian_path_exists = True
    
    # Initialize count for semi-balanced nodes.
    semibalanced = 0
    
    for node in nodes:
        # It is possible for a node not to be a key in the graph dictionary,
        # if it only has incoming edges.
        # In that case: len(graph.get(node, "")) evaluates to 0.
        outdegree = len(graph.get(node, ""))
        indegree = len([value for value in graph.values() if node in value])
        
        if outdegree != indegree: # this node is  not balanced; check if it's semibalanced
            if abs(outdegree - indegree) == 1:
                semibalanced += 1 # this node is semibalanced, keep track
            else: # this node is neither balanced nor semibalanced, so the graph can't contain an eulerian path
                eulerian_path_exists = False
                break
            
            if semibalanced > 2:
                eulerian_path_exists = False
                break
                
    return eulerian_path_exists
    
    

def make_eulerian_graph(graph):
    """
    Add an edge to graph to make it eulerian.
    
    graph: a directed graph that has an eulerian path,
    but not an eulerian cycle. It has two semibalanced nodes,
    and the rest are balanced.
    
    keys: nodes; strings or integers
    values: neighbor nodes; lists of strings or lists of integers
    
    The function finds the start and end node, and then adds
    an edge directed from the end node to the start node.
    
    The outputs are:
    1) the modified graph with the added edge.
    2)  the start node of the original graph
    3) the end node of the original graph.
    """
    if is_eulerian(graph):
        # no conversion needed
        return graph
        
    if not has_eulerian_path(graph):
        print("This graph does not contain an eulerian path.")
        return None
        
    nodes = find_all_nodes(graph)
    
    start = None
    end = None
    
    for node in nodes:
        if node not in graph.keys():
            end = node
        else:
            incoming = [(neighbor, node) for neighbor in graph.keys() if node in graph[neighbor]]
            if len(incoming) ==  0:
                start = node
                
    # Add new edge directed from end to start.
    graph[end] = [start]
    
    # Ensure the modified graph is eulerian
    assert is_eulerian(graph), "Your graph is still not eulerian..."
                
    
    return graph, start, end
    
def find_eulerian_cycle(graph, start = None, used_edges = set()):
    """
    Find the eulerian cycle of a graph, assuming one exists.
    
    graph: directed graph; dictionary where:
    keys: nodes; strings or integers
    values: connected nodes; lists of strings or lists of integers
    
    start: any node with untransversed edges; None or str or int
    used_edges: edges that have been used; set
    
    There are two outputs:
    1) used_edges: all the edges that have been used; list of int or str
    Length should be equal to the total number of edges in the graph.
    2) path: the eulerian cycle; list of integers or strings.
    
    start starts as None, and used_edges starts as an empty set.
    The function firsts finds a balanced subgraph, starting at any node.
    
    If the balanced subgraph is an eulerian cycle, it is returned.
    The used_edges array is also returned.

    Otherwise, the function finds a node that still has unvisited edges
    in the path and sets it to be 'start'. The used_edges array is updated,
    the paths are merged, and the function keeps calling itself with 
    updated values until the result is an eulerian cycle.
    
    """
    
    # Get the total number of edges in the graph; needed for comparison.
    total_edges = sum([len(value) for value in graph.values()])
    
    # Find the first balanced subgraph.
    # Pick start vertex at random.
    current_vertex = None
    
    if start == None:
        current_vertex = choice(list(graph.keys()))
    else:
        current_vertex = start
        
    # Initialize the path.
    path = [current_vertex]
    
    walk = True
    
    # Find the first balanced subgraph.
    while walk:
        neighbors = graph.get(current_vertex)
        if neighbors == None: # this should never happen if all nodes are balanced.
            print("Warning: Your graph is not Eulerian.")
            return [], []
        else:
            neighbor_edges = [(current_vertex, neighbor) for neighbor in neighbors]
            
            # Get the neighbor edges that have not already been used.
            possible_directions = list(set(neighbor_edges) - used_edges)
            
            # If there are no possible directions, you've reached the end of the walk.
            if len(possible_directions) == 0:
                walk = False
            
            # Otherwise, chose a direction and add it to used_edges.
            else:
                edge = choice(possible_directions)
                next_vertex = edge[1] 
                used_edges.add(edge)
                path.append(next_vertex)
                
                current_vertex = next_vertex
                
    # If all edges have been used, the balanced subgraph found is 
    # the eulerian cycle; return it.
    if len(used_edges) == total_edges:
        return used_edges, path
        
    # If not all edges have been used, it means there is still a node w
    # in the path that has untransversed edges.
    
    else:
        # Find w.
        w = None
        for node in path:
            outgoing = [(node, neighbor) for neighbor in graph.get(node) if (node, neighbor) not in used_edges]
            # The incoming thing was mostly a sanity check during development.
            # If the graph is eulerian and there are untransversed outgoing edges, 
            # there should be an equal number of untransversed incoming edges.
            incoming = [(neighbor, node) for neighbor in graph.keys() if (node in graph[neighbor]) and 
            ((neighbor, node) not in used_edges)]
            
            if len(outgoing) > 0 and len(outgoing) == len(incoming):
                w = node
            
                # Set w as the starting point and the updated set of used_edges, and 
                # have the function call itself recursively until it finds the full path.
                
        
    
                new_used_edges, new_path = find_eulerian_cycle(graph, start = w, used_edges = used_edges)
                
    
               # Merge the paths and update the used_edges array
                w_idx = path.index(w)
                path = path[: w_idx] + new_path + path[w_idx + 1:]
                used_edges.update(new_used_edges)
        
    return used_edges, path
    
def make_graph_from_spectrum(spectrum, alphabet):
    """
    Convert a spectrum to a graph.
    
    spectrum: list of strings; each string is a kmer
    alphabet: symbols used to construct the kmers; list of strings
    
    k = length of kmer
    l = k - 1
    
    The graph is a dictionary where
    keys: nodes (l-mers) from where the edge starts; strings
    values: l-mers connected to the key node; lists of strings
    
    A directed edge from lmer1 and lmer2 exists if the suffix
    of lmer1 is the prefix of lmer2.
    """
    graph = {}
    
    k = len(spectrum[0])
    l = k - 1
    
    lmers = ["".join(permutation) for permutation in list(permutations(alphabet, l))]
    
    for lmer in lmers:
        for kmer in spectrum:
            prefix = kmer[: -1]
            suffix = kmer[1: ]
            
            if lmer == prefix:
                if lmer not in graph:
                    graph[lmer] = [suffix]
                else:
                    if suffix not in graph[lmer]:
                        graph[lmer].append(suffix)
            elif lmer == suffix:
                if prefix not in graph:
                    graph[prefix] = [lmer]
                else:
                    if lmer not in graph[prefix]:
                        graph[prefix].append(lmer)
                    
    return graph    
                
        
def find_eulerian_path(graph):
    """
    Return the eulerian path of graph. assuming one exists.
    
    graph: directed graph; dictionary where:
    keys: nodes (strings)
    values: neighbor nodes; lists of strings
    
    start: the semibalanced node at the beginning of the graph
    end: the semibalanced node at the end of the graph
    """
    
    # Make sure graph has an eulerian path
    if not has_eulerian_path(graph):
        return None
        
    # Add an edge between the semibalanced nodes (direction: end -> start)
    # so the graph becomes eulerian
    modified_graph, start, end = make_eulerian_graph(graph)
    
    # Find eulerian cycle in the modified graph
    used_edges, cycle = find_eulerian_cycle(modified_graph)
    
    # Cut the edge.
    start_idx = cycle.index(start)
    end_idx = cycle.index(end)
    
    path = cycle[start_idx: ] + cycle[1: end_idx + 1]
    
    return path
    
def print_eulerian_path(path):
    """
    Print en eulerian path with arrows between the steps.
    
    path: eulerian path; list of strings or integers
    
    An eulerian cycle is a subtype of eulerian path,
    so the function also works for cycles.
    """
    formatted_path = []
    
    last = len(path) - 1
    
    for i, item in enumerate(path):
        formatted_path.append(str(item))
        if i != last:
            formatted_path.append("->")
            
    print(" ".join(formatted_path))
                

if __name__ == "__main__":

    # GRAPH FROM FIG 8.22
    graph_822 = {'J':['D'],'D':['C'],'I':['H'],'H':['F'],'F':['G','E'],\
        'C':['I','A'],'G':['J'],'E':['A'],'A':['F','B'],'B':['C']}
    # A SLIGHTLY BIGGER GRAPH
    bigger_graph = {1:[2], 2:[3], 7:[3,11],\
        4:[5,10],5:[6],6:[3],3:[8,9,1],\
        8:[4],9:[7,4],\
        10:[9],11:[12],12:[7]}
    # SPECTRUM FROM FIG 8.20
    s = ['ATG','TGG','TGC','GTG','GGC','GCA','GCG','CGT']
    
    # ALPHABETS
    bases = ["A", "C", "G", "T"]

    
    # Check is graph_822 is Eulerian.
    
    print("Examining graph_822...")
    if is_eulerian(graph_822):
        print("Graph 822 is eulerian.")
    print("\n")
        
    # Check if graph_822 has an Eulerian path.
    
    print("Examining graph_822...")
    if has_eulerian_path(graph_822):
        print("Graph 822 has an Eulerian path.")
    print("\n")
    
    # Find Eulerian cycle in graph_822
    
    print("Searching for an Eulerian cycle in graph 822...")
    used_edges_822, cycle_822 = find_eulerian_cycle(graph_822)
    print_eulerian_path(cycle_822)
    print("\n")
    
    # Find multiple Eulerian cycles in graph 822
    
    print("Searching for Eulerian cycles in graph 822...")
    print("\n")
    for i in range(3):
        print("Eulerian cycle no: " + str(i + 1))
        used_edges_822, cycle_822 = find_eulerian_cycle(graph_822, None, set())
        print_eulerian_path(cycle_822)
        print("\n")
        
    # Construct graph from spectrum
    
    print("Constructing graph_820 from spectrum s...")
    print("\n")
    graph_820 = make_graph_from_spectrum(s, bases)
    for k, v in graph_820.items():
        print(k, v)
    print("\n")
    
    # Check whether graph 820 1) is Eulerian, and 2) has an Eulerian path.
    
    print("Examining graph 820...")
    print("\n")
    if not is_eulerian(graph_820):
        print("Graph 820 is not eulerian.")
    print("\n")
    if has_eulerian_path(graph_820):
        print("Graph 820 has an eulerian path.")
    print("\n")
    
    # Search for Eulerian path in graph 820
    
    print("Searching for Eulerian path in graph 820...")
    print("\n")
    path_820 = find_eulerian_path(graph_820)
    print("This is the Eulerian path found in graph 820:")
    print("\n")
    print_eulerian_path(path_820)
    print("\n")
    print("This is the sequence it corresponds to:")
    print("\n")
    sequence = []
    last = len(path_820)
    for i, item in enumerate(path_820):
        if i == last - 1:
            sequence.append(item)
        else:
            sequence.append(item[0])
    print("".join(sequence))
    print("\n")
    
    # Find Eulerian cycles in the bigger graph.
   
    print("First, let's check if bigger_graph is eulerian")
    print("\n")
    if is_eulerian(bigger_graph):
        print("The bigger graph is Eulerian.")
    print("\n")
    print("Since the graph is Eulerian, let's try to find an Eulerian cycle.")
    print("\n")
    used_edges_bigger, cycle_bigger = find_eulerian_cycle(bigger_graph, None, set())
    print("This is the eulerian cycle found for bigger_graph: ")
    print("\n")
    print_eulerian_path(cycle_bigger)
    print("\n")
    print("Let's find four different Eulerian cycles in bigger graph:")
    print("\n")
    for i in range(4):
        print("Eulerian cycle no: " + str(i+1))
        print("\n")
        used_edges_bigger, cycle_bigger = find_eulerian_cycle(bigger_graph, None, set())
        print_eulerian_path(cycle_bigger)
        print("\n")
    
        
    
    

    
        
    
        
    
        
    
    
       
    
        
