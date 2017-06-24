# STACS

As a framework for embodiment and external feedback in simulated neural networks, STACS is composed of two main components. First is the neural network [simulator](https://github.com/fywang/stacs/docs/simulator.md) developed using the [Charm++ parallel programming framework](http://charmplusplus.org/). Second is the [interface](https://github.com/fywang/stacs/docs/interface.md) of this simulator to external devices using the [YARP protocol](http://www.yarp.it/). Although these two components interact, neural network simulation may proceed standalone when no external feedback is needed.

## Charm++

An understanding of how programs written in Charm++ is helpful in understanding the programmatic flow of STACS.
Charm++ exists as an extension to the C++ programming language, additionally providing an adaptive runtime system for parallel execution.
The programming model is composed of parallel objects, called _chares_, where computation proceeds according to a message driven paradigm.
That is, communication and the flow of computation in Charm++ occurs by way of asynchronous messages evoking methods on parallel objects.
Here, data dependencies of the evoked methods are satisfied by the communicated messages.
Subsequently, any data generated by these methods are communicated to dependent methods.
Although this message driven paradigm may be achieved in any number of message passing libraries, Charm++ provides an expressive high level framework for developing parallel applications in this manner.

## Persistence

A snapshot of the neural network is composed of details such as network's graph structure, the state of the neural models, and any messages in transit.
This information is stored as files on disk and serializes everything that's needed for the neural network to persist between simulation runs.
An initial snapshot of the neural network is generated by a companion program called [genet](https://github.com/fywang/genet/), whereas STACS is only responsible for reading in a snapshot of the network prior to simulation and writing out a similar snapshot post simulation.

The files that compose the snapshot can be differentiated between network configuration, data files, and records, and are identified by an extension following a base file name for the network.


### Configuration

* `filebase.model` - contains a list of the network models and parameters specific to the particular neural network being simulated
* `filebase.dist` - contains a rolling sum of the amount of data entries per partition of the network (these are the number of vertices, edges, states, and events)


### Data

Because data files contain information that are generally split across multiple computing units, they are further identified by a number corresponding to the computing unit that data resides on.

* `filebase.coord.#` - three dimensional coordinate locations of the vertices in the neural network
* `filebase.adjcy.#` - connection information of the vertices following a condensed sparse row format
* `filebase.state.#` - state information that is relevant to the vertices and their corresponding incoming edges
* `filebase.event.#` - event information that is relevant to the vertices and their corresponding incoming edges

### Records

Records contain information that accumulates over simulation time and are identified by an iteration count spaced evenly by the specified record interval. They are also split across multiple computing units

* `filebase.evtlog.#.#` - event logging
* `filebase.record.#.#` - periodic records

### Groups

These files are used in the persistence and analysis of polychronous groups that emerge from self-organization of the network, and are identified by the *mother* vertex of that group.

* `filebase.group.#` - spatio-temporal stamp of a polychronous group
* `filebase.chart.#` - the above, but with additional information between points

## Execution

STACS is executed by the Charm++ runtime system.

An example execution call is: `./charmrun +p2 ./stacs config.yml`

* `./charmrun` - launches the Charm++ runtime system on the parallel computing platform
* `+p2` - corresponds to launching the Charm++ runtime system with 2 computing units
* `./stacs` - is the program (STACS) to be executed by the Charm++ runtime system
* `config.yml` - is the main configuration file read by STACS

Some key items in the main configuration file are:

* `npdat` - the number of data files there are (this should match with the +p argument for the number of computing units provided to Charm++)
* `npnet` - the number of parallel partitions of the network there are (this should be greater than or equal to npdat)
* `filebase` - the base file name of the neural network snapshot files