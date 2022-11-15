# Setting up an SNN

## Overview

Note: the assumption in the following example(s) will be that the necessary model files are already available and all that's needed are the YAML file configs for defining the network.

There are three YAML files that need to be composed:

- The main config file sets up information related to how the network will be processed by the simulator
- The network graph file defines the structure of the vertices and edges that compose the network
- The network model file defines the substrate dynamics used by the vertices and edges

## Main configuration file

These are files in the `configs` directory. As shown in the example files, there are four main sets of config variables that get set up: simulation, network, timing, and groups.

#### Simulation

The simulation config variables concern with how STACS behaves.

Note: The development of STACS was driven by work on polychronizing networks in a robotics setting, so a number of the available features may not be as relevant for general purpose SNN simulation.

- `runmode`, which is the main relevant variable, can be chosen from the set [simulate, findgroup, and estimate], with the default being "simulate". 
- `randseed` sets a seed for repeatable execution (as long as network partitions are the same).
- `plastic` and `episodic` toggle some additional functionality to the simulation control loop.
- `rpcport` provides the YARP-based port address that STACS will listen to for RPC commands, and `rpcpause` toggles whether the simulator waits for the start command. These variables are ignored if the simulator is not compiled with YARP.

#### Network

The network variables are related to how the network is stored and processed.

- `netwkdir` points to the directory where the data is stored, and `recordir` and `groupdir` are subdirectories in `netwkdir` where records and (polychronous) group information is stored
- `filebase` provides the base name of the network files, and `fileload` and `filesave` are suffixes to `filebase` for which files correspond to which network snapshot to read/write from/to during loading/saving, respectively
- `netfiles` declares how many different files the network is split between on disk, and `netparts` declares how many different partitions the network is split between in simulation

#### Timing

The timing variables are related to how STACS manages the simulation time stepping. All units for the time variables are in milliseconds.

- `tstep` determines how long a simulation time step is, and `tmax` determines the max simulation time
- `teventq` (in conjunction with `tstep`) determines the size of the event queue for incoming events on a vertex
- `tdisplay` sets how often STACS gives an update (printed to STDOUT) on which iteration it's on
- `trecord` determines how often probed network data (e.g. generated events) are recorded to file
- `tsave` determines how often network state (i.e. full snapshot) is saved to disk

#### Groups

These variables are probably not relevant for general simulations, but they are related to working with polychronous groups. There is probably a better way of doing this.

- `grpactive` defines the vertex populations where polychronous groups can be found, and `grpmother` and `grpanchor` provide unique identifiers of a group based off of a certain method of finding them
- `grpminlen` defines the minimum number of vertices in a group, and `grpmaxdur` defines the maximum duration to search over per identifier
- `grpvtxmin` and `grpvtxmax` were used to partition the vertices to find groups over, so that different ranges could be processed in parallel, depending on the resources available.

## Network graph file

These are found in the `networks` directory under a subdirectory related to that specific network (e.g. `networks/polynet`) and should contain the correct `filebase` as specified by the main config file. The network graph file defines the structure of how the different vertices and edges are related to each other. There are also streams for communication into/outof the STACS simulated SNN. 

#### Vertices

This section contains a list of each type of vertex model used by the network, which are uniquely identified by `modname` (this model name should correspond to the model definitions in the network model file).

For any given vertex model, there are some additional instantiation parameters with respect to the graph:

- `order` determines how many of the vertex model to instantiate, thus creating a population of said vertex model. Each of these instances are initialized based on information in the associated network model file.
- `shape` and related parameters such as `radius` and `coord` are used to specify the spatial location of the vertex (or vertices if order > 1).
  - These are useful for spatially dependent connectivity and partitioning. Options for `shape` include "point", "circle", and "sphere"
  - If location information is unnecessary, a "point" at coordinate [0.0, 0.0, 0.0] is the default.

#### Streams

These are essentially a special type of vertex that contains information for real-time I/O through YARP. They are unique (order = 1), so only require a `modname` and `coord` to specify where they are in the network. They are connected to other vertices in the graph in the same way (through edge models).

#### Edges

This section contains a list of each type of edge connecting the different vertex populations. These are uniquely identified by `source` and `target` populations, and how they are connected through an edge model (identified by `modname`).

Note: these define a directed graph, and given `source` and `target` pairing can currently only be specified to be connected with a single edge model.

For any given connection between vertices, we have the following parameters:

- `source` provides the `modname` of the sending vertex model/population (see the above), and `target` provides a list of `modname` of the receiving vertex model(s)/population(s)
- `modname` in the edge specification determines the edge model that is used to compute the connection dynamics (i.e. send/receive processing).
- `cutoff` is used when spatially dependent connectivity is used, and the default is 0.0 for no cutoff distance between any pair of vertices (where beyond that, the `connect` information is not evaluated).
- `connect` contains the information needed to determine whether to create a connection when evaluating a pair of source/target vertices.
	- `conntype` specifies the connection type based on probability distributions, which can currently be one of [uniform, sigmoid, index].
	- Each connection type takes on certain parameters, such as `prob` and `slope`, that specify the probability distribution, or `srcoff` and `srcmul`, that specify index-based offsets/multipliers w.r.t. the source index.

## Network model file

These are also found in the `networks` directory under a subdirectory related to that specific network (e.g. `networks/polynet`) and should contain the correct `filebase` as specified by the main config file. The network model file defines model parameters and state initialization details for the network substrate models (located in the `models` directory). These include model definitions for both vertex and edge models, and they are split up into different YAML "files" that are separated within the file using the `---` and `...` syntax.

The model specifications follow the same general format regardless of whether or not the model pertains to a vertex or an edge (or a stream).

- `type` determines whether the model is for a [vertex, edge, stream]
- `modname` provides the unique model name that is used to differentiate between the different models that compose a given network
- `modtype` provides a reference to the underlying model dynamics that the model uses. Different models (with different `modname`) can reference to the same `modtype` (e.g. this is useful if you have the same type of neuron model that is used in multiple populations).
- `param` specifies model-specific parameters (as implemented by the `modtype` reference) that influence the model dynamics. These parameters are shared across all instances of `modname`, and are defined with of `name` and `value` pairs.
- `state` specifies instantiation details for model-specific state (as implemented by the `modtype` reference) that are either unique per model instantiation and/or are mutable throughout network simulation.
	- These are defined using a `name` as well as the `init` type for initialization:
		- `constant` takes a `value` that all instantiations are initialized with.
		- `uniform` and `normal` and its variants (e.g. over an interval) draw values from a specified probability distribution 
		- `linear` is proportional to distance
	- The `rep` (short for representation) of a given state selects between [real, tick] and determines whether the state is a real number or an integer-based value (used for timestamps)
- `port` specifications are unique to streams, and are used to open communication ports for I/O
