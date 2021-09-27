import numpy as np

class Graph:
    """ A class that constructs a graph of binary skeleton image

    Arguments:
        bin_img: binary image.
    """
    neighbour_grid = np.array([
            [-1,  0],
            [-1, -1],
            [0,  -1],
            [1,  -1],
            [1,   0],
            [1,   1],
            [0,   1],
            [-1,  1]])
    
    class Node:
        """ Simple node clas for storing purposes

        Arguments:
            loc: tuple of row and colum index
        """
        def __init__(self, loc):
            self.r, self.c = loc
            self.neighbours = []
            self.key = str(self.r)+','+ str(self.c)
            self.type = "Node"

    def __init__(self, bin_img):
        self.nodes = {}
        self.endpoints = []
        self.bin_img = bin_img
        self.max_r, self.max_c = bin_img.shape
        self.endpoint_paths = None
        self.amount_clusters = None
        self.__graphify()
        self.__determine_amount_clusters()

    def get_node(self, loc):
        """ Retrieves a node for a given location

        Arguments:
            loc: tuple lof row and column index

        Returns:
            The node
        """
        r, c = loc
        key = str(r) +','+ str(c)

        if key not in self.nodes:
            self.nodes[key] = self.Node(loc) 
        
        return self.nodes[key]

    def __add_neighbours(self, node):
        """ Add neighbours to a node based on the binary image

        Arguments:
            node: The node
        """

        # copy location template
        neighbour_indx = np.copy(Graph.neighbour_grid)

        # apply node location to template neighbourhood
        neighbour_indx[:, 0] += node.r 
        neighbour_indx[:, 1] += node.c

        # Get valid indices and add all neighbours to that node
        valid_indx = self.__get_in_bounds_indx(neighbour_indx)
        for i in range(valid_indx.shape[0]):

            # Check whether valid indx in the binary image is tagged
            r, c = valid_indx[i,:]
            if self.bin_img[r, c] != 0:

                # Get node from dict and add to neighbourhood
                neighbour = self.get_node(valid_indx[i,:])
                node.neighbours.append(neighbour)

        # If node only has one neighbour, tag it as an endpoint
        if len(node.neighbours) < 2:
            node.type = "Endpoint"
            self.endpoints.append(node)

    def __get_in_bounds_indx(self, neighbour_indx):
        """ Checks whether potential indices are in bounds

        Arguments:
            neighbours_indx: 2D array of indices

        Returns:
            1D array of bools (whether in bounds or not)
        """
        r_bool = (neighbour_indx[:, 0] < self.max_r) * (neighbour_indx[:, 0] >= 0)
        c_bool = (neighbour_indx[:, 1] < self.max_c) * (neighbour_indx[:, 1] >= 0)

        return neighbour_indx[(r_bool * c_bool)]


    def __graphify(self):
        """ Helper function to create a graph from the binary image
        """
        # Find all non-zero values in image
        non_zeros = np.where(self.bin_img != 0) # (y, x)
        amount = non_zeros[0].shape[0]

        # For each non-zero values get the node and add neighbours
        for i in range(amount):
            loc = np.array([non_zeros[0][i], non_zeros[1][i]])
            node = self.get_node(loc)
            self.__add_neighbours(node)

        return self.nodes
    
    def __list_intersection(self, list1, list2):
        """ Does an intersection operation between two lists

        Arguments:
            list1: a list
            list2 a list

        Returns:
            As list with all intersecting values in both lists
        """
        temp = set(list2)
        return [val for val in list1 if val in temp]

    def __list_substract(self, list1, list2):
        """ Removes all elements from list1 that are not in list 2

        Arguments:
            list1: a list
            list2 a list

        Returns:
            The new list
        """
        temp = set(list2)
        return [val for val in list1 if val not in temp]
    
    def __remove_node_from_list(self, list, node):

        for i in range(len(list)):
            list_node = list[i]
            if list_node.r == node.r and list_node.c == node.c:
                del list[i]
                break

        return list

    def __determine_amount_clusters(self):
        """ Helper function that determines the number of non-connected subgraphs
        """
        endpoints = self.endpoints.copy()
        # del endpoints[0]
        if len(endpoints) == 0:
            self.amount_clusters = 0
            self.endpoint_paths =  []
            self.max_cluster_length = 0
            return

        cluster_max_lengths = np.zeros(len(endpoints))
        cluster_paths = []
        
        i = 0
        while len(endpoints) > 0:
            
            endpoint = endpoints[0]
            self.__remove_node_from_list(endpoints, endpoint)

            paths = self.get_bfs_shortest_paths(endpoint)
            cluster_paths.append(paths)

            for path in paths:
                other_endpoint = path[-1]
                path_len = len(path)

                if path_len > cluster_max_lengths[i]:
                   cluster_max_lengths[i] = path_len 

                self.__remove_node_from_list(endpoints, other_endpoint)

            i += 1

        self.amount_clusters = i
        self.endpoint_paths = cluster_paths[np.argmax(cluster_max_lengths)]
        self.max_cluster_length = np.max(cluster_max_lengths)

    def get_bfs_shortest_paths(self, starting_endpoint=None):
        """ Retrieve all shortest path from an end point to any other end point

        Arguments:
            starting_point: optional, a node from with all shortest paths should be retrieved,
                otherwise the first one is chosen
        Returns:
            A list of nodes
        """
        # Choose first endpoint as starting point, so subsampling is conistent

        if starting_endpoint == None:
            starting_endpoint = self.endpoints[0]

        explored = []
        queue = [[starting_endpoint]]
        endpoint_paths = []
    
        while queue:
            path = queue.pop(0)
            node = path[-1]

            if node not in explored:
                explored.append(node)
                neighbours = node.neighbours

                for neighbour in neighbours:

                    if neighbour not in explored:
                        new_path = list(path)
                        new_path.append(neighbour)
                        queue.append(new_path)

                        if neighbour.type == "Endpoint":
                            endpoint_paths.append(new_path)
                
        # self.endpoint_paths = endpoint_paths

        return endpoint_paths

    def get_longest_merged_path(self, paths):
        """ Retrieves the longest path in the graph.

        Arguments:
            paths: A list of paths, where each path is a list of nodes

        Returns:
            The longest path, a list of nodes
        """
        longest_path = None
        longest_len = 0

        # find initial longest path
        paths_set = []
        for i in range(len(paths)):
            p_s = paths[i]
            paths_set.append(p_s)
            if len(p_s) > longest_len:
                longest_path = p_s
                longest_len = len(p_s)

        # neighours combindations
        indx_combs = [
            [[0, 0], [0, 0]],
            [[0, 0], [-1, -1]],
            [[-1, -1], [0, 0]],
            [[-1, -1], [-1, -1]]
        ]

        # Search for longest path combination in the given paths
        # Yes, this should be refactored
        for i in range(len(paths)):
            for j in range(len(paths)):
                if i is not j:

                    # Merge paths 
                    intersection = self.__list_intersection(paths_set[i], paths_set[j])
                    p_i = self.__list_substract(paths_set[i], intersection)
                    p_j = self.__list_substract(paths_set[j], intersection)
                    
                    # Find which endpoints of the paths are neighbours by find min distance all possible combinations
                    distances = []
                    for c in indx_combs:
                        distance = np.sqrt(
                            (p_i[c[0][0]].r - p_j[c[0][1]].r)**2 +
                            (p_i[c[1][0]].c - p_j[c[1][1]].c)**2
                        )
                        distances.append(distance)
                    min_index = np.argmin(distances)
                    min_combo = indx_combs[min_index]

                    # Get common neighour of both paths
                    p_i_neighbours = p_i[min_combo[0][0]].neighbours
                    p_j_neighbours = p_j[min_combo[1][0]].neighbours
                    neighbours = self.__list_intersection(p_j_neighbours, p_i_neighbours)
                    try:
                        neighbour = neighbours.pop()
                    except: # no neighbour in common, can't be our combination
                        continue

                    # Flip if some of the paths are reversed
                    if min_combo[0][0] == -1:
                        p_i = reversed(p_i)
                    if min_combo[1][1] == 0:
                        p_j = reversed(p_j)
                    
                    # Stich path together
                    merged = list([*p_j, neighbour, *p_i])

                    if len(merged) > longest_len: # we want to discard subsets
                        longest_path = merged
                        longest_len = len(merged)

        return longest_path

    def nodes_to_numpy(self, nodes):
        """ Transforms a list of nodes to row and column numpy vectors
        
        Arguments:
            nodes: a list of nodes
        
        Returns:
            A tuple:
                1. Row indices as a 1D numpy array
                2. Column indices as a 1D numpy array
        """
        r = np.array([node.r for node in nodes])
        c = np.array([node.c for node in nodes])

        return r, c