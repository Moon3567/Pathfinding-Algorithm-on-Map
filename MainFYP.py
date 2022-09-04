import numpy as np
import os,sys
import folium
import xmltodict as xtd
import webbrowser

#Importing OSM file exported from https://www.openstreetmap.org 
with open('map(London).osm', "rb") as osm_fn: #rb means its a read only file in binary format. 
    map_osm = xtd.parse(osm_fn)['osm'] # using the python lib xmltodict which essentially converts XML files to python dictionary. It can also be done with OSM files.


# Applying bounds to each row to obtain min and max latitude and longitude.
# This is done through parasing, where the OSM data extracts the min and max lat and lon accross all geometrics.
# Bounds are essentialy the "rectangular" areas within the map such as buildings or houses and their given in lat/lon min/ max values as seen
y_max = map_osm['bounds']['@maxlat'] # Create tuple y_max to apply bounds on max lat
y_min = map_osm['bounds']['@minlat'] # Create tuple y_min to apply bounds on min lat
x_max = map_osm['bounds']['@maxlon'] # Create tuple y_max to apply bounds on max lon
x_min = map_osm['bounds']['@minlon'] # Create tuple y_min to apply bounds on min lon
parsed_bounds = [x_min, x_max, y_min, y_max] # Store all parsed bounds into "bounds"

######################################################################################################################################################################

# A node is one of the core elements within the OSM file, where it includes a single point of space that is defined by its lat, lon and node id.
# Parsing nodes from OSM file to extract the length of node, lat and log.
Node=map_osm['node'] # Accessing the node key from map_osm and assigning that to "Node"
node_length=len(Node) # Taking the length of that "Node" object.
Nodeid = [0]*node_length # "*" creates duplicates of that list and appends them togerther.

xy = [] # Starts an empty list

# iterates from 0 to the length of the "Node" object
for node_coord in range(node_length): # for example, if that length was 8, "0" through "7" meaning the next value from 0 to 7 is assigned to the "node_coord" variable.
    Nodeid[node_coord]=float(Node[node_coord]['@id']) # returns id for node and turn into floating point number to make it into a simpler dictionary to read
    x=float(Node[node_coord]['@lat']) # Getting lat coordinates of that node and storing into node_coord
    y=float(Node[node_coord]['@lon']) # Getting lon coordinates of that node and storing into node_coord
    xy.append([x,y])
parsed_node={'id':Nodeid, 'xy':xy} # Applying to one dict with "id" key pointing to a list of node ID'S and the "xy" key pointing to lat-long coordinates

######################################################################################################################################################################
#"Way" is one of the fundamentals elements of the map which is an ordered list of nodes. Essentially it is the path or rather a line that represents a linear feature such as roads, river, wall. 
# Within this section, it looks into the OSM dictionary and checks to see if certain values are there, converts them into appropriate pyton type and selects specific fields to be put back into a dictionary
way_map=map_osm['way'] # Accessing the "way" key from map_osm and assigning it to "Paths"

Node_ways=len(way_map) # Taking the length of "way"

Wayid=[0]*Node_ways # "*" creates duplicates of that list 
nodes_path=[0]*Node_ways # "*" creates duplicates of that list 
tags=[0]*Node_ways # "*" creates duplicates of that list 

# In this section the same will be done as above with bounds.
for way in range(Node_ways):
  temp_way = way_map[way] # get key id "way" from OSM data and assign it to "temp_way" variable
  Wayid[way] = float(temp_way['@id']) # return all "way" id into floating point number to make it into a simpler dictionary to read
  # 'nd' is a feature which is constructed by references to the node id.
  reference_node=len(temp_way['nd']) # gets the length of the temp_way and assigns it to reference node
  temp_reference_node=[0]*reference_node # "*" creates duplicates of that list to st
  
  # essentially gets "temp_reference_node", inserts it into list and checks for another referenced node, to be returned using the float() method. 
  for reference in range(reference_node):
    temp_reference_node[reference]=float(temp_way['nd'][reference]['@ref']) # referencing another node and returning it as float point number so it can be put into the dictionary
  
  nodes_path[way] = temp_reference_node # assigning the "temp_reference_node" variable into node path variable which is the "way" tag. This will be the fully parsed variable for node paths

# Essentially a loop which searches for certain key-value pairs so it can be converted into a python dictionary.
# The .key() method is used to view objects that contains keys of the dictionary, as a list.
  if 'tag' in temp_way.keys(): # The key id 'tag' defines what the "way" represent. For example if its a motorway, road, etc.
    if type(temp_way['tag']) is list: # if the 'tag' for that specific "way" is within the list 
      tags[way]=temp_way['tag'] # add to list and this will be sent to be returned to a floating point number 
    else:
          tags[way]=[temp_way['tag']] # same as above
  else:
    tags[way]=[] # empty list to store "tag"

parsed_way={'id':Wayid,'nodes':nodes_path, 'tags':tags} # 'id' represents the wayid, 'node' represents the node path and 'tags' represents the "way" keys and values


######################################################################################################################################################################
# The relation is an element of ordered list. This list is include multiple nodes, ways and/or relations. For example several roads for bus routes
# The same process has been done to parse the "relation" key id extracted from the OSM file
relation_map=map_osm['relation'] # Accessing the node key from map_osm and assigning that to "relation"
node_relation=len(relation_map) # length of the relation between the node or the way
Relation_id=[0]*node_relation # "*" creates duplicates of that list 

for relation in range(node_relation): 
    currentRelation = relation_map[relation] #  The current relation between a node or way
    currentId=currentRelation['@id'] # Gets the id of that current relation
    Relation_id[relation]=float(currentId)  # converts currentid into floating point number and assigns to "Relation_id" variable. This help put the relation id into a simpler dictionary to read
    
parsed_relation={'id':Relation_id} # Using "id" key which points to a list of "relation_id"

######################################################################################################################################################################

#Assigning all parsed variables into "parsed_osm"
parsed_osm={ 
    'bounds':parsed_bounds, # uses bound variable to fetch the bounds key from OSM file
    'way':parsed_way, # uses parsed_way variable to fetch the way key from OSM file
    'node':parsed_node, # uses parsed_node variable to fetch the node key from OSM file
    'relation':parsed_relation, # uses parsed_relation variable to fetch the relation key from OSM file
    'attributes':map_osm.keys() # fetches the attributes from OSM file  
}

ways_num = len(parsed_way['id']) # getting length of the way id from OSM file
ways_node_set = parsed_way['nodes'] # setting each node to an id
node_ids = dict() #empty dict
length_nodeid = len(parsed_node['id']) # getting length of the node id

for i in range(length_nodeid):
    node_ids[parsed_node['id'][i]] = i # storing i inside the empty dict which essentially stores the id of that node
    

#https://wiki.openstreetmap.org/wiki/Key:highway#Roads
# The OSM file contains many types of roads such as, moterways, residential ways, highway, etc. Essentially selecting the type of road to be displayed on the map
road_vals = [
        "highway",
        "motorway",
        "motorway_link",
        "trunk",
        "trunk_link",
        "primary",
        "primary_link",
        "secondary",
        "secondary_link",
        "tertiary",
        "road",
        "living_street",
        "motorway_junction",
    ]
                
######################################################################################################################################################################

#Building an Adjacency matrix of vertex to vertex (node to node) allows it to become adjacent or rather a neighbour if there is a path to between both nodes.  
def create_matrix():
    adjacency_matrix = np.full((node_length,node_length), float('inf')) # "node_length" used twice for vertex to vertex and make the matrix infinate
    np.fill_diagonal(adjacency_matrix, 0) # fill the main diagonal of the given array of any dimensionality and fill with 0's
    
    # for loop is used to say if the tags are included within the variable road_val then dont skip 
    for currentWay in range(ways_num): #returns the length of the way id
        skip = True 
        for tags in parsed_way['tags'][currentWay]: # tags in key id "way"
            if tags['@k'] in road_vals: # if the tags are included within the variable "road_vals" then dont skip
                skip = False
                break
        if skip: # if skip is false
            continue
        node_set = ways_node_set[currentWay] # set node on map 
        nodes_num=len(node_set) # length of the set node

# In this section the adjacency matrix checks the boolean value to see if there is a direct path between two verticies or in this case two nodes. 

        for firstnode_local_index in range(nodes_num): # Checks the first node by the node id
            firstnode_id = node_set[firstnode_local_index] # The variable "firstnode_id" is made to allow the system to check both the id of the node and length
            firstnode_index = node_ids.get(firstnode_id, -1) # Get id for node
            if firstnode_index == -1: continue # If the value = -1 then continue
            
#  The same has been done as above but this time it checks to see if the new node that has been selected, can make a path from the first node to the second.

            for othernode_local_index in range(firstnode_local_index + 1, nodes_num): 
                othernode_id=node_set[othernode_local_index]
                othernode_index = node_ids.get(othernode_id, -1)
                if othernode_index==-1: continue # If the value = -1 then continue

# Essentially in order excecute this code there must be two conditions to be met. the first being the variables firstnode_id MUST equal othernode_id
# The second being the value from adjacency_matrix MUST be infinate. If either of those conditions are not met, the code will not excecute. 
                if(firstnode_id != othernode_id and adjacency_matrix[firstnode_index,othernode_index]==float('inf')):
                    adjacency_matrix[firstnode_index, othernode_index] = 1 # If above meets conditions, make value equal to 1
                    adjacency_matrix[othernode_index, firstnode_index] = 1 # If above meets conditions, make value equal to 1

    return adjacency_matrix

######################################################################################################################################################################

# In this section, the Dijkstra's Algorithm is implemented using the adjacency matrix (vertex to vertex) to find the shortest path
# The ith element in "paths" is the next step along the shortest path starting at vertex "node"
def dijkstra(source_node, adjacency_matrix, paths): # "source_node" is the starting node
    seen = dict() # The "seen" dictionary keeps track of which vertices have been visited so far.
    seen[source_node] = True # the function set "seen[source_node]" is set to true
    paths[source_node] = source_node # storing the source nodes 

    n_nodes = len(adjacency_matrix)# length of adjacency matrix
    current_node = source_node # This is the target node
    current_distance = float('inf') # infinate as the adjacency matrix has no direct edge between the node and neighbour node
    
    # Within this for loop, the algorithm finds the closest immediate neighbour to the source node
    for node in range(n_nodes): # for a node that is in range within the adjacency matrix
        node_distance = adjacency_matrix[source_node][node] # check distance between that node against the source node
        # by going through all the nodes that have not been visited
        if node != source_node and node_distance < current_distance: # if that node is not equal to the source node and the node distance is less than current distance
            current_node = node # make that node the current node
            current_distance = node_distance # make the node distance, the current distance

    
    # Within the while loop, the algorithm loops through every node within the adjacency matrix to find the neighbours of the current node
    while node > 0: # distance from source node from itself is always 0
        previous_node  = source_node
        temp_distance = float('inf')
    
        for neighbours in range(n_nodes):# for each unvisited neighbour node within the adjacency matrix
            # the algorithm compares the old cost of current value of the shortest path from the start vertex to the neighbour
            if seen.get(neighbours, False) == False and adjacency_matrix[source_node][current_node] != float('inf') and adjacency_matrix[current_node][neighbours] != float('inf'):
                
                # the variable "cost" updates the "seen" dict if there is an edge from source node to neighbour node
                cost = adjacency_matrix[source_node][current_node] + adjacency_matrix[current_node][neighbours] 
                adjacency_matrix[source_node][neighbours] = min(adjacency_matrix[source_node][neighbours], cost) # choose the vertex with minimum distance from source to neighbour 
                adjacency_matrix[neighbours][source_node] = adjacency_matrix[source_node][neighbours] #  distance from neighbour to source 
            
            # In this section, if the new shortest path has been found, then updates list and says THIS is the shortest path
                if adjacency_matrix[source_node][neighbours] == cost: # the value of the shortest path from the distance between the neighbour and current node is the new cost
                    paths[neighbours] = current_node
                elif adjacency_matrix[source_node][neighbours] == 1: # new cost is 1, where the value of the shortest path from the distance between the neighbour and current node
                    paths[neighbours] = source_node
                    
                # if the source node and neighbour node is less than the old cost (temp_distance), put the neighbour and its cost to the priority queue and update the list where the shortest path is kept.
                if adjacency_matrix[source_node][neighbours] < temp_distance: 
                    previous_node = neighbours
                    current_distance = adjacency_matrix[source_node][neighbours]

        if previous_node == source_node: break # if the shortest path has been found
        seen[previous_node] = True # finished with this node
        current_node = previous_node
        node -= 1 # subtract right operand from left and then assign to left operand. If true then both are equal
        
######################################################################################################################################################################

# Plotting the routes onto the interactive map by creating a continious node path.
def plot_routes(seen, adjacency_matrix):
    paths = dict() # empty dictionary
    dijkstra(seen, adjacency_matrix, paths) # using the dijkstra's function to store the adjacency_matrix and the 2 dictionaries

    node_values = [] # empty list that will store the plotted nodes
    for plot_node in paths.keys(): # .key() method used which essentially reflect any changes in paths
        add_node = [plot_node, 0] # the "add_node" will add the node paths into the list and since its 0, it means it will start with 0 nodes added
        while paths[plot_node] != plot_node: # while the list is empty
            add_node[1] += 1 # Add one node to the list
            plot_node = paths[plot_node] 
        node_values.append(add_node) 
        

    return node_values, paths
######################################################################################################################################################################

# In this section there will be 3 maps to be generated, the first map will display all nodes within the map.
# The second map will display the closest nodes connected to the primary node selected by user
# The final map will display a route from the source node to the destination node 

def PlotAllNodes():
    x0, y0 = (float(parsed_bounds[2]), float(parsed_bounds[0])) # Variable x0 and y0 represent the x-axis and corresponding y-axis of bounds within the map
    x1, y1 = (float(parsed_bounds[3]), float(parsed_bounds[1])) # Same as above
    center = ((x0+x1)/2, (y0+y1)/2) # calculates the center of which the nodes are displayed
    first_map = folium.Map(location = center, zoom_start = 20) # positions the map in the center of where nodes are plotted

    for i in range(length_nodeid):
        # fetching nodes from the variable "parsed_node" and returning an iterable by assigning the index 0 to x and 1 to y
        xy = (parsed_node['xy'][i][0], parsed_node['xy'][i][1]) #getting the coordinates of all nodes by using OSM id 'xy'
        folium.CircleMarker(xy, radius=3, color="green", fill=True, fill_color="green", popup=str(i)).add_to(first_map) # Set all nodes to colour green 
    return first_map

######################################################################################################################################################################
# Second Map 
def Closest_Nodes(source_node, node_values):
    x0, y0 = (float(parsed_bounds[2]), float(parsed_bounds[0])) # Variable x0 and y0 represent the x-axis and corresponding y-axis of bounds within the map
    x1, y1 = (float(parsed_bounds[3]), float(parsed_bounds[1])) #  Same as above
    
    center = ((x0+x1)/2, (y0+y1)/2) # calculates the center of which the nodes are displayed
    second_map = folium.Map(location = center, zoom_start = 20) # positions the map in the center of where nodes are plotted

    for i,j in node_values: # fetching a list of all plotted nodes
        
        xy = (parsed_node['xy'][i][0], parsed_node['xy'][i][1]) #getting the coordinates of all nodes by using OSM id 'xy'
        
        if(i != source_node): # If the plotted nodes is NOT equal to the source node 
            
            folium.CircleMarker(xy, radius=3, color="purple", fill=True, fill_color="green", popup=str(i)).add_to(second_map) # make that node red (these will be the closest nodes)
        else:
            folium.CircleMarker(xy, radius=3, color="darkblue", fill=True, fill_color="green", popup=str(i)).add_to(second_map) # or else make it blue (this will be the source node)
    return second_map

######################################################################################################################################################################

#Generating a map to display the path between source and destination


def Destination_Path(i,path):
    
    node_path = [(parsed_node['xy'][i][0], parsed_node['xy'][i][1])] #getting the coordinates of all nodes by using OSM id 'xy' and storing it in "node_path"
    
    while path[i] != i: # while the position value at the ith position in "path" is not equal to "i"
        # add that node into the list (As in if the source node and destination node have been found, then append by getting the xy coorinates of both nodes)
        node_path.append((parsed_node['xy'][path[i]][0], parsed_node['xy'][path[i]][1])) 
        i = path[i] # put both the source node and destination node within the "i" list
        
    third_map = folium.Map(location = node_path[-1], zoom_start = 20) # -1 means the last element in the sequence

    folium.CircleMarker(node_path[-1], radius=6, color="green", fill=True, fill_color="purple").add_to(third_map) # Source node
    folium.Marker(node_path[0], icon = folium.Icon(color="blue", icon="circle", prefix='fa')).add_to(third_map) # Destination node
    
    folium.PolyLine(locations = node_path, weight=3, color="red", dash_array=5).add_to(third_map) # The polyLine to connect the source and destination node
    
    return third_map

######################################################################################################################################################################

# In this section, the three maps will be saved as an .html to view via browser
# the "OpenHTMLMapinBrowser" function is used to allow to open html files via browser 
def OpenHTMLMapinBrowser(filename):
    url = "file://" + os.path.realpath(filename)
    webbrowser.open(url, new=2)

# Call "PlotAllNodes" function to display all nodes on the map
map_1 = PlotAllNodes()
map_1.save("PlotAllNodes.html")
OpenHTMLMapinBrowser("PlotAllNodes.html")


#Third Map Generator to show path from source to destination
while(True):
    SourceNode = int(input("Enter Node ID:")) # Only accepts int values
    adjacency_matrix = create_matrix() # Assigns the function "create_matrix" to "adjacency_matrix" variable
    # the "plot_routes" functions is used to plot the source node and all nodes which are within the adjacency_matrix to "node_values" and "p"
    node_values, path = plot_routes(SourceNode, adjacency_matrix)
    
    if(not SourceNode): # type 0 to end
        print("Ended")
        sys.exit(1)

    
    map_2 = Closest_Nodes(SourceNode, node_values) # Display source node and nodes the closest nodes
    map_2.save("Closest_Nodes.html")
    OpenHTMLMapinBrowser("Closest_Nodes.html")

    while(True):
        DestinationNode = int(input("Enter the DestinationNode ID (0=exit, -1=new route): ")) # After entering Node ID then input destination node ID
        
        if DestinationNode == -1: # Enter -1 to exit the map
            break
        
        if(not DestinationNode):
            print("Ended") # Print after -1 has been entered
            sys.exit(1)
            
        map_3 = Destination_Path(DestinationNode, path) # Display the destination node and the source node.
        map_3.save("Destination_Path.html")
        OpenHTMLMapinBrowser("Destination_Path.html")