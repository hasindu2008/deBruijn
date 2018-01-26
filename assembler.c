/*
Cython module used to interface to assembler routines.
*/
#include <stdlib.h>
#iunclude <stdio.h>
#include <string.h>
#include <math.h>
#include "common.c>"

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

int REF = 1;
int READ = 2;
int REF_AND_READ = 3;

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// Encapsulate information related to a node of the de-Bruijn graph,
// i.e. a sequence k-mer and related data.
//
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// Forward declaration of Edge struct, so it can be referred to in Node.
struct Edge;
struct EdgeStack;

// Represents a node in the graph
struct Node:
    Edge* edges[4];
    char* sequence;
    int colours;
    int position;
    int kmerSize;
    int nEdges;
    double weight;
    char dfsColour;

// Simple implementation of a stack, for storing nodes.
struct NodeStack:
    Node** elements;
    int capacity;
    int top;

// Represents an edge in the graph
struct Edge:
    Node* startNode;
    Node* endNode;
    double weight;

// A stack of edges
struct EdgeStack:
    Edge** elements;
    int capacity;
    int top;

// A dictionary of Nodes
struct NodeDict:
    Node*** buckets;
    int* bucketSize;
    int nBuckets;

// Hold a path through the graph
struct Path:
    NodeStack* nodes;
    int nNodes;
    int isBubble;
    double weight;

// A stack of paths
struct PathStack:
    Path** elements;
    int capacity;
    int top;

// A graph
struct DeBruijnGraph:
    int kmerSize;
    NodeStack* allNodes;
    NodeDict* nodes;

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

Path* createPath(int initialSize):
    /*
    Create and return a Path struct.
    */
    Path* thePath = <Path*>(malloc(sizeof(Path)));

    if thePath == NULL:
        logger.error("Could not allocate path");

    thePath.nodes = createNodeStack(initialSize);
    thePath.nNodes = 0;
    thePath.isBubble = 0;
    thePath.weight = 0.0;

    return thePath;

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void destroyPath(Path* thePath):
    /*
    free up memory in Path struct.
    */
    destroyNodeStack(thePath.nodes);
    thePath.nodes = NULL;
    free(thePath);
    thePath = NULL;

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void addNodeToPath(Path* thePath, Node* theNode, double weight):
    /*
    Add a Node to the specified path.
    */
    if thePath == NULL:
        logger.error("Null path");

    if theNode == NULL:
        logger.error("Null Node");

    NodeStack_Push(thePath.nodes, theNode);
    thePath.nNodes += 1;
    thePath.weight += weight;

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

Node* createNode(char* sequence, int colour, int position, int kmerSize, double weight):
    /*
    Create and return a new node.
    */
    Node* theNode = <Node*>(malloc(sizeof(Node)));

    if theNode == NULL:
        logger.error("Could not allocate node");

    // Need to allocate kmerSize + 1 to store the null terminating character
    theNode.edges[0] = NULL;
    theNode.edges[1] = NULL;
    theNode.edges[2] = NULL;
    theNode.edges[3] = NULL;

    theNode.sequence = sequence;
    theNode.weight = weight;
    theNode.colours = colour;
    theNode.dfsColour = 'N';
    theNode.kmerSize = kmerSize;
    theNode.position = position;
    theNode.nEdges = 0;

    return theNode;

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void* destroyNode(Node* theNode):
    /*
    free all memory associated with a Node struct. The Node takes ownership
    of all associated string data.
    */
    // Edges belong to the node, so destroy these now.
    int i = 0;

    for i in range(theNode.nEdges):
        destroyEdge(theNode.edges[i]);

    free(theNode);
    theNode = NULL;

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

inline int Node_Equal(Node* thisNode, Node* otherNode):
    /*
    Two nodes are equal if and only if their sequencese are equal.
    */
    if thisNode.kmerSize != otherNode.kmerSize:
        return 0;
    else:
        return (strncmp(thisNode.sequence, otherNode.sequence, thisNode.kmerSize) == 0);

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

inline void Node_AddEdge(Node* theNode, Edge* theEdge):
    /*
    Add the specified edge to the specified node.
    */
    theNode.edges[theNode.nEdges] = theEdge;
    theNode.nEdges += 1;

    if theNode.nEdges > 4:
        logger.error("Node struct cannot have > 4 edges. Now we have %s. This will cause problems" %(theNode.nEdges));

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


inline char* Node_GetSuffix(Node* theNode):
    /*
    Return the suffix string of this node. The suffix is simply the last kmerSize-1
    elements of the node sequence.
    */
    return theNode.sequence + 1;

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


inline char* Node_GetPrefix(Node* theNode):
    /*
    Return the suffix string of this node. The suffix is simply the first kmerSize-1
    elements of the node sequence.
    */
    return theNode.sequence;

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


inline int nodePosComp(const void* x, const void* y):
    /*
    Comparison function for Node structs, for use in qsort, to sort Nodes by their
    positions.
    */
    const Node** nodeOne = <const Node**>(x);
    const Node** nodeTwo = <const Node**>(y);

    // Sort by position;
    return  nodeOne[0].position - nodeTwo[0].position;

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


inline unsigned int hashKmer(char* theKmer, int size):
    /*
    Return a hash value for the specified kmer.
    */
    int i = 0;
    unsigned int hashVal = 0;

    // 1) Sum integers formed by 4-character sub-strings.
    for i in range(0, size-4, 4):
        hashVal += theKmer[i];
        hashVal += (theKmer[i + 1] << 8);
        hashVal += (theKmer[i + 2] << 16);
        hashVal += (theKmer[i + 3] << 24);

    // 2) multiply by 101 and add new value. From Paul Larson (see StackOverflow)
    //for i in range(size):
    //    hashVal = hashVal * 101 + theKmer[i]

    return hashVal;

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

PathStack* createPathStack(int capacity):
    /*
    Create and return a stack for storing Paths.
    */
    PathStack* theStack = <PathStack*>(malloc(sizeof(PathStack)));

    if theStack == NULL:
        logger.error("Error. Could not allocate path stack with capacity %s" %(capacity));

    theStack.elements = <Path**>(malloc(sizeof(Path*)*capacity));

    if theStack.elements == NULL:
        logger.error("Error. Could not allocate path stack elements with capacity %s" %(capacity));

    theStack.top = -1 // Empty;
    theStack.capacity = capacity;
    return theStack;

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void destroyPathStack(PathStack* theStack):
    /*
    Clears up memory in stack. Does not destroy the
    Paths.
    */
    int i = 0;

    for i in range(theStack.top + 1):
        destroyPath(theStack.elements[i]);

    free(theStack.elements);
    theStack.elements = NULL;
    free(theStack);
    theStack = NULL;

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


int PathStack_IsEmpty(PathStack* theStack):
    /*
    Return 1 if the stack is empty and 0 otherwise.
    */
    if theStack.top == -1:
        return 1;
    else:
        return 0;

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


int PathStack_IsFull(PathStack* theStack):
    /*
    Return 1 if the stack is full and 0 otherwise.
    */
    if theStack.top + 1 == theStack.capacity:
        return 1;
    else:
        return 0;

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void PathStack_Push(PathStack* theStack, Path* element):
    /*
    Add a new element to the stack. Elements always go on the top,
    i.e. in the highest position. Realloc if necessary.
    */
    Path** temp = NULL

    // Need to realloc
    if PathStack_IsFull(theStack):
        temp = <Path**>(realloc(theStack.elements, 2*sizeof(Path*)*theStack.capacity))

        if temp == NULL:
            logger.error("Could not re-allocate PathStack");
        else:
            theStack.elements = temp;
            theStack.capacity *= 2;

    theStack.top += 1;
    theStack.elements[theStack.top] = element;

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

Path* PathStack_Pop(PathStack* theStack):
    /*
    Pop and return the top element from the stack.
    */
    Path* thePath = NULL

    if PathStack_IsEmpty(theStack):
        return NULL;
    else:
        thePath = theStack.elements[theStack.top]
        theStack.top -= 1;
        return thePath;

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

NodeStack* createNodeStack(int capacity):
    /*
    Create and return a stack for storing nodes.
    */
    NodeStack* theStack = <NodeStack*>(malloc(sizeof(NodeStack)));

    if theStack == NULL:
        logger.error("Error. Could not allocate node stack with capacity %s" %(capacity));

    theStack.elements = <Node**>(malloc(sizeof(Node*)*capacity));

    if theStack.elements == NULL:
        logger.error("Error. Could not allocate node stack elements with capacity %s" %(capacity));

    theStack.top = -1 // Empty;
    theStack.capacity = capacity;
    return theStack;

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void destroyNodeStack(NodeStack* theStack):
    /*
    Clears up memory in stack. Does not destroy the
    nodes.
    */
    free(theStack.elements);
    theStack.elements = NULL;
    free(theStack);
    theStack = NULL;

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


int NodeStack_IsEmpty(NodeStack* theStack):
    /*
    Return 1 if the stack is empty and 0 otherwise.
    */
    if theStack.top == -1:
        return 1;
    else:
        return 0;

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


int NodeStack_IsFull(NodeStack* theStack):
    /*
    Return 1 if the stack is full and 0 otherwise.
    */
    if theStack.top + 1 == theStack.capacity:
        return 1;
    else:
        return 0;

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void NodeStack_Push(NodeStack* theStack, Node* element):
    /*
    Add a new element to the stack. Elements always go on the top,
    i.e. in the highest position. Realloc if necessary.
    */
    Node** temp = NULL;

    // Need to realloc
    if NodeStack_IsFull(theStack):
        temp = <Node**>(realloc(theStack.elements, 2*sizeof(Node*)*theStack.capacity));

        if temp == NULL:
            logger.error("Could not re-allocate NodeStack");
        else:
            theStack.elements = temp;
            theStack.capacity *= 2;

    theStack.top += 1;
    theStack.elements[theStack.top] = element;

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

Node* NodeStack_Pop(NodeStack* theStack):
    /*
    Pop and return the top element from the stack.
    */
    Node* theNode = NULL

    if NodeStack_IsEmpty(theStack):
        return NULL;
    else:
        theNode = theStack.elements[theStack.top]
        theStack.top -= 1;
        return theNode;

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// TODO: Write function for printing nodes.
//    def __repr__(self):
//        /*
//        Returns a string representation of the node.
//        */
//        return "%s (%s). Pos = %s" %(self.sequence, self.colours, self.position)
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// Structs and functions for handling Edges.
// An Edge struct encapsulates an edge of the de-Bruijn graph.
//
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

Edge* createEdge(Node* startNode, Node* endNode, double weight):
    /*
    Create and return a new edge.
    */
    Edge* theEdge = <Edge*>(malloc(sizeof(Edge)));

    if theEdge == NULL:
        logger.error("Error. Could not allocate edge");

    theEdge.startNode = startNode;
    theEdge.endNode = endNode;
    theEdge.weight = weight;

    return theEdge;

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void destroyEdge(Edge* theEdge):
    /*
    free memory used by Edge struct.
    */
    free(theEdge);
    theEdge = NULL;

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

EdgeStack* createEdgeStack(int capacity):
    /*
    Create and return a stack for storing Edges.
    */
    EdgeStack* theStack = <EdgeStack*>(malloc(sizeof(EdgeStack)));

    if theStack == NULL:
        logger.error("Error. Could not allocate edge stack with capacity %s" %(capacity));

    theStack.elements = <Edge**>(malloc(sizeof(Edge*)*capacity));

    if theStack.elements == NULL:
        logger.error("Error. Could not allocate edge stack elements with capacity %s" %(capacity));

    theStack.top = -1 // Empty
    theStack.capacity = capacity;
    return theStack;

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void destroyEdgeStack(EdgeStack* theStack):
    /*
    Clears up memory in stack. Does not destroy the
    Edges.
    */
    free(theStack.elements);
    theStack.elements = NULL;
    free(theStack);
    theStack = NULL;

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


int EdgeStack_IsEmpty(EdgeStack* theStack):
    /*
    Return 1 if the stack is empty and 0 otherwise.
    */
    if theStack.top == -1:
        return 1;
    else:
        return 0;

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


int EdgeStack_IsFull(EdgeStack* theStack):
    /*
    Return 1 if the stack is full and 0 otherwise.
    */
    if theStack.top == theStack.capacity - 1:
        return 1;
    else:
        return 0;

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void EdgeStack_Push(EdgeStack* theStack, Edge* element):
    /*
    Add a new element to the stack. Elements always go on the top,
    i.e. in the highest position. Realloc if necessary.
    */
    Edge** temp = NULL;

    // Need to realloc
    if EdgeStack_IsFull(theStack):
        temp = <Edge**>(realloc(theStack.elements, 2*sizeof(Edge*)*theStack.capacity));

        if temp == NULL:
            raise StandardError, "Could not re-allocate EdgeStack";
        else:
            theStack.elements = temp;
            theStack.capacity *= 2;

    theStack.top += 1;
    theStack.elements[theStack.top] = element;

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

Edge* EdgeStack_Pop(EdgeStack* theStack):
    /*
    Pop and return the top element from the stack.
    */
    Edge* theEdge = NULL;

    if EdgeStack_IsEmpty(theStack):
        return NULL;
    else:
        theEdge = theStack.elements[theStack.top];
        theStack.top -= 1;
        return theEdge;

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

NodeDict* createNodeDict(int nBuckets):
    /*
    Create and return a dictionary of kmer/void pointers.
    */
    NodeDict* theDict = <NodeDict*>(malloc(sizeof(NodeDict)));

    if theDict == NULL:
        logger.error("Error. Could not NodeDict");

    theDict.buckets = <Node***>(malloc(nBuckets*sizeof(Node**)));

    if theDict.buckets == NULL:;
        logger.error("Error. Could not NodeDict.buckets of size %s" %(nBuckets));

    theDict.bucketSize = <int*>(malloc(nBuckets*sizeof(int)));

    if theDict.bucketSize == NULL:
        logger.error("Error. Could not NodeDict.bucketSize of size %s" %(nBuckets));

    theDict.nBuckets = nBuckets;

    int i = 0;

    for i in range(nBuckets):
        theDict.buckets[i] = NULL;
        theDict.bucketSize[i] = 0;

    return theDict;

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void destroyNodeDict(NodeDict* theDict):
    /*
    free memory used by NodeDict.
    */
    int i = 0;
    int j = 0;

    for i in range(theDict.nBuckets):
        free(theDict.buckets[i]);

    free(theDict.buckets);
    free(theDict.bucketSize);
    free(theDict);
    theDict = NULL;

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

int NodeDict_FindOrInsert(NodeDict* theDict, Node** theNode, int keyLen, Node** nodeForUpdating):
    /*
    Either return the element which is associated with the specified key, or
    insert a new element at the relevant position. Return 0 if the element was not found,
    and 1 if it was.
    */
    int hashValue = hashKmer(theNode[0].sequence, keyLen) % theDict.nBuckets;
    int bucketSize = 0;
    int i = 0;
    int initialBucketSize = 5;
    Node* testNode = NULL;
    Node* newNode = NULL;

    // Need to allocate new bucket
    if theDict.buckets[hashValue] == NULL:
        theDict.buckets[hashValue] = <Node**>(malloc(initialBucketSize*sizeof(Node*)));

        if theDict.buckets[hashValue] == NULL:
            logger.error("Could not allocate hash table bucket with size %s" %(initialBucketSize));

        // Always set to NULL first
        for i in range(initialBucketSize):
            theDict.buckets[hashValue][i] = NULL;

        newNode = createNode(theNode[0].sequence, theNode[0].colours, theNode[0].position, theNode[0].kmerSize, theNode[0].weight);
        theDict.buckets[hashValue][0] = newNode;
        theDict.bucketSize[hashValue] = initialBucketSize;
        theNode[0] = newNode;
        return 0;

    // Bucket is there. Check all elements in bucket for match.
    else:
        bucketSize = theDict.bucketSize[hashValue];

        for i in range(bucketSize):

            // Found empty slot. Insert new element
            if theDict.buckets[hashValue][i] == NULL:
                newNode = createNode(theNode[0].sequence, theNode[0].colours, theNode[0].position, theNode[0].kmerSize, theNode[0].weight);
                theDict.buckets[hashValue][i] = newNode;
                theNode[0] = newNode;
                return 0;

            // Check for match
            else:
                // Match. Return this element.
                testNode = theDict.buckets[hashValue][i];
                if theNode[0] == testNode or strncmp(theNode[0].sequence, testNode.sequence, keyLen) == 0:
                    nodeForUpdating[0] = testNode;
                    return 1;

    // If we get here, then we found no empty slots, and no matches. So we need to
    // allocate more space in the relevant bucket and then insert the key/value pair
    // in the next available space..
    int oldBucketSize = theDict.bucketSize[hashValue];
    int newBucketSize = 2*oldBucketSize;
    Node** temp = <Node**>(realloc(theDict.buckets[hashValue], sizeof(Node*)*newBucketSize));

    if temp == NULL:
        raise StandardError, "Could not re-allocate bucket";
    else:
        // Set new entries to NULL
        for i in range(oldBucketSize, newBucketSize):
            temp[i] = NULL;

        newNode = createNode(theNode[0].sequence, theNode[0].colours, theNode[0].position, theNode[0].kmerSize, theNode[0].weight);
        theDict.bucketSize[hashValue] = newBucketSize;
        theDict.buckets[hashValue] = temp;
        theDict.buckets[hashValue][oldBucketSize] = newNode;
        theNode[0] = newNode;

    return 0;

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

DeBruijnGraph* createDeBruijnGraph(int kmerSize, int nBuckets):
    /*
    Allocate memory for graph data.
    */
    DeBruijnGraph* theGraph = <DeBruijnGraph*>(malloc(sizeof(DeBruijnGraph)));

    if theGraph == NULL:
        logger.error("Could not allocate memory for DeBruijnGraph");

    theGraph.kmerSize = kmerSize;
    theGraph.allNodes = createNodeStack(nBuckets);
    theGraph.nodes = createNodeDict(nBuckets);

    return theGraph;

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void destroyDeBruijnGraph(DeBruijnGraph* theGraph):
    /*
    Free memory used by graph.
    */
    int i = 0;

    // Destroy all nodes
    for i in range(theGraph.allNodes.top + 1):
        destroyNode(theGraph.allNodes.elements[i]);

    // These only hold pointers to memory which will be
    // freed elsewhere.
    destroyNodeStack(theGraph.allNodes);
    destroyNodeDict(theGraph.nodes);
    free(theGraph);

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

int DeBruijnGraph_InsertOrUpdateNode(DeBruijnGraph* theGraph, Node** theNode):
    /*
    Check if a node is already present. If it is, update it, otherwise
    insert it. Return 1 if the node was inserted and 0 if it was updated.
    */
    Node* nodeForUpdating = NULL;
    int foundNode = NodeDict_FindOrInsert(theGraph.nodes, theNode, theGraph.kmerSize, &nodeForUpdating);

    // Need to create a new node, copying values from theNode
    if not foundNode:
        NodeStack_Push(theGraph.allNodes, theNode[0]);
        return 1;

    // Update existing node
    else:
        // Update colours of nodes already in graph.
        nodeForUpdating.colours |= theNode[0].colours;
        nodeForUpdating.weight += theNode[0].weight;
        theNode[0] = nodeForUpdating;
        return 0;

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void DeBruijnGraph_AddEdge(DeBruijnGraph* theGraph, Node* startNode, Node* endNode, double weight):
    /*
    */
    double oldWeight = startNode.weight;
    int startNodeWasInserted = DeBruijnGraph_InsertOrUpdateNode(theGraph, &startNode);
    int endNodeWasInserted = DeBruijnGraph_InsertOrUpdateNode(theGraph, &endNode);
    Edge* newEdge = NULL;

    int i = 0;

    // Check all outgoing edges from startNode. If it has none, then make one for this edge. Otherwise,
    // check all existing edges for a match, and update accordingly.
    for i in range(4):
        if startNode.edges[i] == NULL:
            newEdge = createEdge(startNode, endNode, weight);
            Node_AddEdge(startNode, newEdge);
            break;
        elif startNode.edges[i].endNode == endNode:
            startNode.edges[i].weight += weight;
            break;
        else:
            continue;
    else:
        pass; // This only happens when there are Ns in the sequence.
        //logger.error("Error in assembler. Could not find matching end-node for edge. Something is very wrong")
        //logger.error("Start node sequence is %s. End node sequence is %s" %(startNode[0].sequence[0:startNode[0].kmerSize], endNode[0].sequence[0:endNode[0].kmerSize]))
        //raise StandardError, "Assembly Error!!"

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

int dfsVisit(Node* theNode, double minWeight):
    /*
    */
    int nEdges = theNode.nEdges;
    int i = 0;
    Node* nextNode = NULL;
    Edge* edge = NULL;

    theNode.dfsColour = 'g';

    for i in range(nEdges):
        edge = theNode.edges[i];

        // Ignore low-weight edges that are only in reads
        if edge.endNode.colours == READ and edge.weight < minWeight:
            continue;

        nextNode = edge.endNode;

        if nextNode.dfsColour == 'w':

            // Found cycle in this path
            if dfsVisit(nextNode, minWeight) == 1:
                return 1;
            // This path ok. Go to next edge
            else:
                continue;
        elif nextNode.dfsColour == 'g':
            // Found cycle
            //logger.debug("Found cycle. From %s to %s" %(theNode.position, nextNode.position))
            return 1;

        // Black node. Already explored past this node. Go to next edge.
        else:
            continue;

    // No cycles in any path reachable from this node.
    theNode.dfsColour = 'b';
    return 0;

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

int detectCyclesInGraph_Recursive(DeBruijnGraph* theGraph, double minWeight):
    /*
    */
    Node** allNodes = theGraph.allNodes.elements;
    Node* thisNode = NULL;
    Node* nextNode = NULL;
    int i = 0;
    int j = 0;
    int nNodes = theGraph.allNodes.top + 1;
    int foundCycle = 0;

    for i in range(nNodes):
        thisNode = allNodes[i];
        thisNode.dfsColour = 'w';

    for i in range(nNodes):
        thisNode = allNodes[i];

        if thisNode.dfsColour == 'w':
            // Found cycle
            foundCycle = dfsVisit(thisNode, minWeight);

            if foundCycle == 1:
                return True;

    return False;

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

int detectCyclesInGraph(DeBruijnGraph* theGraph, double minWeight):
    /*
    */
    Node* thisNode;
    Node* nextNode;
    Edge* edge;
    int i = 0;
    int j = 0;
    int nEdges = 0;
    int nNodes = theGraph.allNodes.top + 1;
    Node** allNodes = theGraph.allNodes.elements;

    int nFilledBuckets = 0;
    int nEntriesThisBucket = 0;

    for i in range(theGraph.nodes.nBuckets):
        if theGraph.nodes.buckets[i] != NULL:
            nFilledBuckets += 1;
            for j in range(theGraph.nodes.bucketSize[i]):
                if theGraph.nodes.buckets[i][j] != NULL:
                    nEntriesThisBucket += 1;
                else:
                    break;

    logger.debug("nNodes = %s. nFilledBuckets = %s. mean entries/bucket = %s" %(nNodes, nFilledBuckets, float(nEntriesThisBucket)/nFilledBuckets));
    qsort(<void*>allNodes, nNodes, sizeof(Node*), nodePosComp);

    for i in range(nNodes):
        thisNode = allNodes[i];
        thisNode.dfsColour = 'w';

    Node* sourceNode = allNodes[0];
    Node* endNode = allNodes[nNodes-1];
    NodeStack* theStack = createNodeStack(nNodes);
    int reachedEnd = False;

    NodeStack_Push(theStack, sourceNode);

    while not NodeStack_IsEmpty(theStack):

        thisNode = NodeStack_Pop(theStack);

        if thisNode.dfsColour == 'w':
            thisNode.dfsColour = 'g';
        elif thisNode.dfsColour == 'g':
            thisNode.dfsColour = 'b';
        else:
            pass;


        nEdges = thisNode.nEdges;

        for i in range(nEdges):
            edge = thisNode.edges[i];
            nextNode = edge.endNode;

            // TODO: temp hack. Replace with Nodes_Equal later
            if Node_Equal(nextNode, endNode):
                reachedEnd = True;

            if nextNode.dfsColour == 'w':
                NodeStack_Push(theStack, nextNode);

            // Found a cycle
            elif nextNode.dfsColour == 'g':
                destroyNodeStack(theStack);
                return True;

            else:
                pass;

    destroyNodeStack(theStack);
    return False;

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

char* createSequenceFromPath(Path* thePath):
    /*
    Create and return a string representation of the sequence of a specific path
    through the graph.
    */
    int nNodes = thePath.nNodes;
    char* theString = <char*>(malloc( (nNodes+1)*sizeof(char)));

    if theString == NULL:
        logger.error("Could not allocate memory for string of size %s" %(nNodes+1));

    int i = 0;

    for i in range(nNodes):
        theString[i] = thePath.nodes.elements[i].sequence[0];

    theString[nNodes] = 0;
    return theString;

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

int checkPathForCycles(Path* thePath):
    /*
    Check if this path contains a cycle.
    */
    int nNodes = thePath.nNodes;
    int i = 0;

    //logger.debug("Checking path with %s nodes for cycles" %(nNodes))

    // Set all dfs colours to white
    for i in range(nNodes):
        thePath.nodes.elements[i].dfsColour = 'w';

    // Check all nodes in order. If we see the same node twice, then
    // there is a cycle. If we get to the end without seeing any nodes twice,
    // then no cycle.
    for i in range(nNodes):
        if thePath.nodes.elements[i].dfsColour == 'w':
            thePath.nodes.elements[i].dfsColour = 'g';
        else:
            //logger.debug("Found cycle")
            return 1;

    //logger.debug("No cycles")
    return 0;

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

PathStack* getVariantPathsThroughGraphFromNode(DeBruijnGraph* theGraph, Path* thePath, double minWeight):
    /*
    Check all valid paths through the graph starting at the last node in "thePath". If any path
    returns to the reference sequence, then stop and add that path to the returned list. Also, if the
    path never returns to the reference, but is sufficiently long, then add this to the list.
    */
    PathStack* thePathStack = createPathStack(10);
    PathStack* finishedPaths = createPathStack(10);
    Path* pathSoFar = NULL;
    Node* endSoFar = NULL;
    Node* newEnd = NULL;
    Edge* theEdge = NULL;
    int nEdgesThisNode = 0;
    int i = 0;
    int j = 0;
    int hasCycle = 0;

    PathStack_Push(thePathStack, thePath);

    while not PathStack_IsEmpty(thePathStack):

        pathSoFar = PathStack_Pop(thePathStack);
        endSoFar = pathSoFar.nodes.elements[pathSoFar.nNodes-1];

        // TODO: Replace with maxHaplotypes??
        if thePathStack.top + 1 > 20 or finishedPaths.top + 1 > 20:
            //logger.info("Too many paths %s (%s) from this node. Giving up" %(thePathStack.top, finishedPaths.top))
            destroyPath(pathSoFar);
            destroyPathStack(thePathStack);
            destroyPathStack(finishedPaths);
            return NULL;

        // At the moment, don't allow any cycles.
        hasCycle = checkPathForCycles(pathSoFar);

        if hasCycle:
            destroyPath(pathSoFar);

        // Got back to ref, this path is done with. This is a bubble.
        elif endSoFar.colours == REF_AND_READ:
            pathSoFar.isBubble = 1;
            PathStack_Push(finishedPaths, pathSoFar);

        // No reads here. Not quite sure how this could happen. Went from ref and read to only
        // ref.
        elif endSoFar.colours == REF:
            destroyPath(pathSoFar);

        // Keep extending path
        else:
            // Dumb check for loops
            //if pathSoFar.nNodes > 1000:

                // If we are at the end of this path, and we haven't returned to the reference yet, then still suggest a variant
                // using all the sequence so far, but only if the path is long enough. This is not a bubble.
                //if nEdgesThisNode == 0:
                //    destroyPath(pathSoFar)
                //    if pathSoFar.nNodes > 50:
                //        pathSoFar.isBubble = 0
                //        PathStack_Push(finishedPaths, pathSoFar)
                //    else:
                //        destroyPath(pathSoFar)

                //else:
            nEdgesThisNode = endSoFar.nEdges;

            for i in range(nEdgesThisNode):
                theEdge = endSoFar.edges[i];
                newEnd = theEdge.endNode;

                if theEdge.weight >= minWeight or newEnd.colours == REF_AND_READ or newEnd.colours == REF:
                    newPath = createPath(theGraph.kmerSize);

                    // Copy old path
                    for j in range(pathSoFar.nNodes):
                        addNodeToPath(newPath, pathSoFar.nodes.elements[j], 0.0);

                    // Weight for this path is weight of existing path + weight of new edge
                    newPath.weight = pathSoFar.weight;
                    addNodeToPath(newPath, newEnd, theEdge.weight);
                    PathStack_Push(thePathStack, newPath);

            destroyPath(pathSoFar);

    destroyPathStack(thePathStack);
    return finishedPaths;

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

list findBubblesInGraph(DeBruijnGraph* theGraph, double minWeight, char* refSeq, char* chrom, int refStart, int refEnd, int assemStart, int assemEnd, int verbosity):
    /*
    Find and return all bubbles in the graph.

    1) Check all nodes
    2) If node has edge going to non-ref node, then follow it
    3) Keep going until we get back to reference.
    4) Keep a stack of the current path
    */
    int nNodes = theGraph.allNodes.top + 1;
    int nEdgesThisNode = 0;
    int nBubblePaths = 0;
    int nNodesThisPath = 0;
    int i = 0;
    int j = 0;
    Edge* theEdge = NULL;
    Node** allNodes = theGraph.allNodes.elements;
    Node* theNode = NULL;
    PathStack* bubblePathsThisNode = NULL;
    Path* thePath = NULL;
    Path* theBubblePath = NULL;
    list variants = [];
    Variant theVar;

    for i in range(nNodes):
        theNode = allNodes[i];

        // First node which is ref only
        if theNode.colours == REF_AND_READ and theNode.position >= assemStart and theNode.position < assemEnd:
            nEdgesThisNode = theNode.nEdges;

            for j in range(nEdgesThisNode):
                //logger.debug("Node %s has %s edges" %(i, nEdgesThisNode))
                theEdge = theNode.edges[j];

                // New kmer in reads. Start path using a depth-first
                // search from this node. Include first node in path.
                if theEdge.endNode.colours == READ:

                    thePath = createPath(theGraph.kmerSize);
                    addNodeToPath(thePath, theNode, 0.0);
                    addNodeToPath(thePath, theEdge.endNode, 0.0);
                    bubblePathsThisNode = getVariantPathsThroughGraphFromNode(theGraph, thePath, minWeight);

                    if bubblePathsThisNode != NULL:
                        nBubblePaths = bubblePathsThisNode.top + 1;

                        //logger.debug("There are %s bubble paths from node %s and edge %s in region %s:%s-%s" %(nBubblePaths, i, j, chrom, assemStart, assemEnd))

                        for j in range(nBubblePaths):
                            theBubblePath = bubblePathsThisNode.elements[j];
                            theVar = extractVarFromBubblePath(theGraph, theBubblePath, refSeq, chrom, refStart, refEnd, verbosity);

                            //if verbosity >= 3:
                            //    logger.debug("Assembler found variant %s. Path has %s nodes, with total weight %s. Weight per node = %s" %(theVar, thePath.nNodes, thePath.weight, thePath.weight/thePath.nNodes))

                            if theVar is not None:
                                variants.append(theVar);

                        destroyPathStack(bubblePathsThisNode);

    return variants;

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void logPath(Path* thePath, char* refSeq, int refStart):
    /*
    Log a path through the graph.
    */
    int i = 0;
    int startPos = thePath.nodes.elements[0].position;
    Node* theNode;
    logger.debug("Logging path of %s nodes" %(thePath.nNodes));

    for i in range(thePath.nNodes):
        theNode = thePath.nodes.elements[i];
        logger.debug("Pos = %s. Seq = %s. RefSeq = %c. Colours = %s. Node weight = %s" %(theNode.position, theNode.sequence[0: theNode.kmerSize], refSeq[startPos + i - refStart], theNode.colours, theNode.weight));

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

Variant extractVarFromBubblePath(DeBruijnGraph* theGraph, Path* bubblePath, char* refSeq, char* chrom, int refStart, int refEnd, int verbosity):
    /*
    Construct and return a representation of the variant implied by this bubble in the graph.
    */
    int nNodesThisPath = bubblePath.nodes.top + 1;
    Node* bubbleStart = bubblePath.nodes.elements[0];
    Node* bubbleEnd = bubblePath.nodes.elements[nNodesThisPath-1];
    Node* theNode;
    int bubbleStartPos = 0;
    int bubbleEndPos = 0;

    // Bubble comes back to reference, so take start and end points on reference 
    // from start and end nodes in bubble
    if bubblePath.isBubble:
        bubbleStartPos = bubbleStart.position;
        bubbleEndPos = bubbleEnd.position;

        if bubbleEndPos < bubbleStartPos:
            //if verbosity >= 3:
                //logger.debug("Found wonky candidate with end pos < start pos. (%s --> %s)" %(bubbleStartPos, bubbleEndPos))
                //logPath(bubblePath, refSeq, refStart)
            bubbleEndPos = bubbleStartPos + bubblePath.nNodes - 1;
            return None;

        //logger.debug("Assembled path with %s nodes. Start = %s. End = %s" %(bubblePath.nNodes, bubbleStart.position, bubbleEnd.position))
        //logPath(bubblePath, refSeq, refStart)

    // Not a bubble, so does not come back to reference. This only happens for long insertions. Use
    // start point in ref as insertion point.
    else:
        bubbleStartPos = bubbleStart.position;
        bubbleEndPos = bubbleStartPos;

    // TODO: We should be able to deal with these now. Just add nNodes to bubble
    if bubbleEndPos < bubbleStartPos:
        logger.warning("Found complex variation that Platypus can't deal with yet at %s:%s-%s" %(chrom, min(bubbleStartPos,bubbleEndPos), max(bubbleStartPos,bubbleEndPos)));
        return None;

    char* charReadSeq = <char*>(malloc(sizeof(char)*(nNodesThisPath+1)));

    if charReadSeq == NULL:
        logger.error("Could not allocate memory for charReadSeq of size %s" %(nNodesThisPath+1));

    int i = 0;

    // Create sequence of variant. Remember to add null terminating element
    for i in range(nNodesThisPath):
        charReadSeq[i] = bubblePath.nodes.elements[i].sequence[0];

    charReadSeq[nNodesThisPath] = 0;

    bytes thisRefSeq = refSeq[bubbleStartPos - refStart: bubbleEndPos - refStart + 1];
    bytes thisReadSeq = charReadSeq;

    // The bytes variable holds a copy, so free the original.
    free(charReadSeq);

    // Push as far to left as possible
    while len(thisReadSeq) > 0 and len(thisRefSeq) > 0:

        if thisRefSeq[-1] == thisReadSeq[-1]:
            thisRefSeq = thisRefSeq[:-1];
            thisReadSeq = thisReadSeq[:-1];
        else:
            break;

    // Trim leading reference sequence
    while len(thisReadSeq) > 0 and len(thisRefSeq) > 0:

        if thisRefSeq[0] == thisReadSeq[0]:
            bubbleStartPos += 1;
            thisRefSeq = thisRefSeq[1:];
            thisReadSeq = thisReadSeq[1:];
        else:
            break;

    //if verbosity >= 3:
    //    logger.debug("After trimming, candidate variant is %s:%s-%s %s --> %s" %(chrom, bubbleStartPos, bubbleEndPos, thisRefSeq, thisReadSeq))

    if verbosity >= 3 and abs(len(thisRefSeq) - len(thisReadSeq)) > 100:
        logger.info("Platypus assembler found candidate variant of size %s. %s:%s-%s %s --> %s" %(len(thisReadSeq) - len(thisRefSeq), chrom, bubbleStartPos, bubbleEndPos, thisRefSeq, thisReadSeq));

    int varLen = len(thisReadSeq) - len(thisRefSeq);
    int varStartPos = bubbleStartPos;

    //if varLen == 0:
    //    thisRefSeq = thisRefSeq[:-1]
    //    thisReadSeq = thisReadSeq[:-1]

    //if verbosity >= 3:
    //    logger.debug("After final final trimming, candidate variant is %s:%s-%s %s --> %s" %(chrom, bubbleStartPos, bubbleEndPos, thisRefSeq, thisReadSeq))

    //if abs(len(thisRefSeq) - len(thisReadSeq)) > 100:
    //    logger.debug("Assembler found large variant %s:%s, %s --> %s. Length = %s" %(chrom, varStartPos, thisRefSeq, thisReadSeq, abs(len(thisRefSeq) - len(thisReadSeq))))

    return Variant(chrom, varStartPos, thisRefSeq, thisReadSeq, 0, ASSEMBLER_VAR);

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void loadReferenceIntoGraph(DeBruijnGraph* theGraph, char* refSeq, int refStart, int kmerSize):
    /*
    Load k-mers from the specified reference sequence into the
    graph.
    */
    int i = 0;
    int lenRef = strlen(refSeq);
    Node tempStartNode;
    Node tempEndNode;

    for i in range( (lenRef-kmerSize) - 1):

        tempStartNode.sequence = refSeq + i;
        tempStartNode.kmerSize = kmerSize;
        tempStartNode.colours = REF;
        tempStartNode.position = refStart + i;
        tempStartNode.weight = 1;

        tempEndNode.sequence = refSeq + i + 1;
        tempEndNode.kmerSize = kmerSize;
        tempEndNode.colours = REF;
        tempEndNode.position = refStart + i + 1;
        tempEndNode.weight = 1;

        DeBruijnGraph_AddEdge(theGraph, &tempStartNode, &tempEndNode, 1);

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

char* createReverseComplementSequence(char* theSeq, int seqLen):
    /*
    Create and return a string containing the reverse complement of the specified
    sequence.
    */
    char* newSeq = <char*>(malloc(sizeof(char)*(seqLen+1)));
    int i = 0;

    for i in range(seqLen):
        if theSeq[i] == 'A':
            newSeq[seqLen - i - 1] = 'T';
        elif theSeq[i] == 'T':
            newSeq[seqLen - i - 1] = 'A';
        elif theSeq[i] == 'C':
            newSeq[seqLen - i - 1] = 'G';
        elif theSeq[i] == 'G':
            newSeq[seqLen - i - 1] = 'C';
        else:
            newSeq[seqLen - i - 1] = theSeq[i];

    newSeq[seqLen] = 0;
    return newSeq;

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void loadReadIntoGraph(cAlignedRead* theRead, DeBruijnGraph* theGraph, int minQual, int kmerSize):
    /*
    */
    char* theSeq = theRead.seq;
    char* theQuals = theRead.qual;
    int length = theRead.rlen;
    int startPos = theRead.pos;
    int i = 0;
    int j = 0;
    int thisMinQual = 100000000;
    int NsInKmer = 0;
    Node tempStartNode;
    Node tempEndNode;

    for i in range( (length-kmerSize) - 1):
        thisMinQual = 100000000;
        NsInKmer = 0;

        for j in range(i, i+kmerSize+1):

            thisMinQual = min(thisMinQual, theQuals[j]);

            if theSeq[j] == 'N':
                NsInKmer = 1;

        if thisMinQual >= minQual and NsInKmer == 0:

            tempStartNode.sequence = theSeq + i;
            tempStartNode.kmerSize = kmerSize;
            tempStartNode.colours = READ;
            tempStartNode.position = - 1;
            tempStartNode.weight = thisMinQual;

            tempEndNode.sequence = theSeq + i + 1;
            tempEndNode.kmerSize = kmerSize;
            tempEndNode.colours = READ;
            tempEndNode.position = -1;
            tempEndNode.weight = thisMinQual;

            DeBruijnGraph_AddEdge(theGraph, &tempStartNode, &tempEndNode, thisMinQual);

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void loadBAMDataIntoGraph(DeBruijnGraph* theGraph, list readBuffers, int assembleBadReads, int assembleBrokenPairs, int minQual, int kmerSize):
    /*
    Load k-mers from the specified BAM file into the graph. K-mers containing
    Ns are ignored, as are k-mers containing low-quality bases.
    */
    bamReadBuffer theBuffer;

    for theBuffer in readBuffers:

        readStart = theBuffer.reads.windowStart;
        readEnd = theBuffer.reads.windowEnd;
        badReadStart = theBuffer.badReads.windowStart;
        badReadEnd = theBuffer.badReads.windowEnd;
        brokenReadStart = theBuffer.brokenMates.windowStart;
        brokenReadEnd = theBuffer.brokenMates.windowEnd;

        while readStart != readEnd:
            if not Read_IsQCFail(readStart[0]):
                loadReadIntoGraph(readStart[0], theGraph, minQual, kmerSize);

            readStart += 1;

        if assembleBadReads:
            while badReadStart != badReadEnd:
                if not Read_IsQCFail(badReadStart[0]):
                    loadReadIntoGraph(badReadStart[0], theGraph, minQual, kmerSize);

                badReadStart += 1;

        if assembleBrokenPairs:
            while brokenReadStart != brokenReadEnd:
                if not Read_IsQCFail(brokenReadStart[0]):
                    loadReadIntoGraph(brokenReadStart[0], theGraph, minQual, kmerSize);

                brokenReadStart += 1;

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

list assembleReadsAndDetectVariants(char* chrom, int assemStart, int assemEnd, int refStart, int refEnd, list readBuffers, char* refSeq, options):
    /*
    */
    int minQual = options.minBaseQual;
    int minMapQual = options.minMapQual;
    int largestVariant = options.assemblyRegionSize;
    int minReads = options.minReads;
    int assembleBadReads = options.assembleBadReads;
    int assembleBrokenPairs = options.assembleBrokenPairs;
    int kmerSize = options.assemblerKmerSize;
    int minWeight = minReads*minQual;
    int nBuckets = 5000;
    int verbosity = options.verbosity;
    list theVars = [];

    if verbosity >= 3:
        logger.debug("Assembling region %s:%s-%s" %(chrom, assemStart, assemEnd));

    DeBruijnGraph* theGraph = createDeBruijnGraph(kmerSize, nBuckets);

    loadReferenceIntoGraph(theGraph, refSeq, refStart, kmerSize);
    loadBAMDataIntoGraph(theGraph, readBuffers, assembleBadReads, assembleBrokenPairs, minQual, kmerSize);

    // If this is true, then don't allow cycles in the graph.
    if options.noCycles:
        while detectCyclesInGraph_Recursive(theGraph, minWeight):
            if kmerSize > 50:
                if verbosity >= 3:
                    logger.debug("Could not assemble region %s:%s-%s without cycles. Max k-mer size tried = %s" %(chrom, assemStart, assemEnd, kmerSize));
                break;
            else:
                if verbosity >= 3:
                    logger.debug("Found cycles in region %s:%s-%s with kmer size %s. Trying again with kmer size %s" %(chrom, assemStart, assemEnd, kmerSize, kmerSize+5));

                kmerSize += 5;
                destroyDeBruijnGraph(theGraph);
                theGraph = createDeBruijnGraph(kmerSize, nBuckets);
                loadReferenceIntoGraph(theGraph, refSeq, refStart, kmerSize);
                loadBAMDataIntoGraph(theGraph, readBuffers, assembleBadReads, assembleBrokenPairs, minQual, kmerSize);
        else:
            theVars = findBubblesInGraph(theGraph, minWeight, refSeq, chrom, refStart, refEnd, assemStart, assemEnd, verbosity);
    else:
        theVars = findBubblesInGraph(theGraph, minWeight, refSeq, chrom, refStart, refEnd, assemStart, assemEnd, verbosity);

    destroyDeBruijnGraph(theGraph);
    //logger.debug("Finished assembling region %s:%s-%s" %(chrom, start, end))
    //logger.debug("Found vars %s" %(theVars))
    return sorted(theVars);

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
