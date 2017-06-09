# Simulator

The simulator proceeds in three main phases: initialization from disk, simulation of the neural network, and finalization to disk.
On disk, snapshots of the neural network enable persistence between simulation runs.
In addition to saving the network state at the end of a simulation run, periodic checkpointing of the network also occurs during simulation.

The following is slightly outdated and needs to be updated to reflect more recent changes to the programmatic flow.

## Initialization

During the initialization phase, STACS processes the configuration files for the neural network and sets up the appropriate parallel objects (chares) across the computing units.
It then reads the snapshot files from disk and starts the simulation after the neural network data is loaded into the appropriate chares.
The relevant function calls and messages that are involved in this phase are as follows.

### Main Entry Point

Programs written in Charm++ have a main entry point defined as the constructor of the "main chare".
For STACS, this is simply named `Main`, and the entry point of the program is at `Main::Main(CkArgMsg *msg)`, where the `CkArgMsg` message provides the standard command line arguments to STACS as `msg->argc` and `msg->argv`.
STACS expects that there is only one command line argument, which is the location (relative path) of the main configuration file.

Reading in the configuration files proceeds as follows:

* `Main::ParseConfig()` - parses the main configuration file, `config.yml`, and sets global readonly variables across the parallel computing platform
* `Main::ReadDist()` - reads in the network configuration file, `filebase.dist`
* `Main::ReadModel()` - reads in the network configuration file, `filebase.model`

The network information is then packaged into messages as:

* `Main::BuildDist()` - builds an `mDist` message, which provides information during construction of the `NetData` chare array
* `Main::BuildModel()` - builds an `mModel` message, which provides information during construction of the `Network` chare array
These messages are then passed to the constructors of the `NetData` and `Network` chare arrays, respectively.

### Constructing Parallel Objects

When an entry method is evoked on the chare array, all chares in the chare array excecute the evoked method.
There are two main groups of parallel objects utilized by STACS, realized as chare arrays in Charm++.
Briefly, the `NetData` chare array is responsible for the network data de/serialization from/to disk, and the `Network` chare array is responsible for the neural network simulation.

During construction of these chare arrays, the main setup procedures involve allocating the relevant data structures according to the configuration of the neural network.
Bookkeeping for which chare of the `NetData` array corresponds to which chare of the `Network` array is also calculated.

#### NetData Chare Array

`NetData::NetData(mDist *msg)` is the constructor for the `NetData` chare array, where the `mDist` message provides information needed to read in the network snapshot.

After the the relevant data structures have been allocated, the network snapshot is read in from disk. Here, `NetData::ReadCSR` reads in the network data files `filebase.coord`, `filebase.adjcy`, `filebase.state`, and `filebase.event`.

#### Network Chare Array

`Network::Network(mModel *msg)` is the constructor for the `Network` chare array, where the `mModel` message provides information regarding the network models and their parameters that will be used during simulation.

### Loading Network Data

After the `NetData` and `Network` chare arrays have been initialized, they return control to `Main::InitSim()`, which acts as a coordination point before the next stage of the initialization phase occurs.
Here, the data structures have been allocated, and data from `NetData` is loaded into `Network`.

The process of moving network data from `NetData` to `Network` proceeds as:

* `Main::InitSim()` - sends a reference to the `NetData` chare array to the `Network` chare array so that it may request network data
* `Network::LoadNetwork(CProxy_NetData cpdat)` - requests its specific network part from `NetData`
* `NetData::LoadNetwork(int prtidx, const CkCallback &cb)` - processes the request from `Network` and sends back the requested network part as a `mPart` message
* `Network::LoadNetwork(mPart *msg)` - receives the network part from `NetData` and populates its data structures accordingly

Once the network data has been loaded into the `Network` chare array, `Network::CreateGroup()` sets up the event communication structure of the network simulation according to the adjacency information.
This communication structure provides multicast sections of the `Network` chare array.

## Simulation

After the network data has been loaded from `NetData` to `Network`, control is returned back to `Main::StartSim()`, which acts as a coordination point before simulation starts.
STACS enters into the simulation phase here.

### Main Simulation Loop

`Network::Cycle()` is the main simulation loop that determines the forward progression of the neural network simulation as well as the periodic export of network data.
There are two main counters that are incremented as the simulation progresses, the iteration count, `iter`, and the simulation time elapsed, `tsim`.

The control flow per cycle checks for:

* `tsim >= tmax` - if the maximum simulation time has been reached
* `iter == checkiter` - if the simulation should checkpoint a network snapshot to disk
* `iter == reciter` - if the simulation should save record data to disk

If the above conditions pass through, the simulation enters into computation of the next iteration or time step of the network, generating event information as well as any simulation records.

This network data is then processed as follows:

* `Network::BuildEvent()` - builds an `mEvent` message, which contains event information to be communicated to other, dependent network parts
* `Network::StoreRecord()` - maintains a local copy of records until they may be written out to disk

At the end of an iteration the `mEvent` messages generated by a network part are communicated to its multicast section of the `Network` chare array.

### Event Communication

`Network::CommEvent(mEvent *msg)` is the main entry method that processes event communication, passed as `mEvent` messages, between network parts.

Once all dependent messages for the next iteration have been received, the network part enters into its next cycle.

### Checking Network Data

STACS periodically saves snapshots of the network data during simulation such that the simulation may be recovered and restarted in case of interruption.

The checkpointing process proceeds as follows:

* `Network::CheckNetwork()` - moves control out of the main simulation loop to perform checkpointing
  * `Network::BuildPart()` - builds an `mPart` message, which contains all the necessary information for serialization of a network partition
  * The network partition information is sent to `NetData`, and control is returned back to the main simulation loop
* `NetData::CheckNetwork(mPart *msg)` - receives network part data from `Network`
 * `NetData::WriteCSR()` - serializes the network data onto disk
* Information about the network distribution is reduced to `Main`
* `Main::CheckNetwork(CkReductionMsg *msg)`, where the msg is of type `netDist`, which contains information about the distribution of the network across the network partitions
 * `Main::WriteDist()` - writes network distribution information to disk

### Checking Network Records

STACS periodically writes records generated from the simulation onto disk so that may be used for offline processing.

The recording process proceeds as follows:

* `Network::CheckRecord()` - moves control out of the main simulation loop to perform checkpointing
  * `Network::BuildRecord()` - builds an `mRecord` message, which contains records information packages record information
  * The records information is sent to the `NetData` chare array, and control is returned back to the main simulation loop
* `NetData::CheckRecord(mRecord *msg)` - receives and processes record data from the `Network` chare array
  * `NetData::WriteRecord()` - serializes the record data onto disk

## Finalization

When the simulation is complete, the `Network` chare array returns control to `Main::StopSim()`, which acts as a coordination point for the network to have finished any computation and communication processes.

### Saving Network Data

The final state of neural network is saved to disk.
This is very similar to checkpointing the network.

The saving process proceeds as follows:

* `Network::SaveNetwork()` - is called to perform the final checkpoint
  * `Network::BuildPart()` - builds an `mPart` message, which contains all the necessary information for serialization of a network partition
  * The network partition information is sent to `NetData`
* `NetData::SaveNetwork(mPart *msg)` - receives network part data from `Network`
  * `NetData::WriteCSR()` - serializes the network data onto disk
* Information about the network distribution is reduced to `Main`
* `Main::SaveNetwork(CkReductionMsg *msg)`, where the msg is of type `netDist`, which contains information about the distribution of the network across the network partitions
  * `Main::WriteDist()` - writes network distribution information to disk

### Saving Network Records

Any records that haven’t been saved yet are written out to disk as well.

The recording process proceeds as follows:

* `Network::SaveRecord()` - is called to perform the final saving of records
  * `Network::BuildRecord()` - builds an `mRecord` message, which contains records information packages record information
  * The records information is sent to the `NetData` chare array
* `NetData::SaveRecord(mRecord *msg)` - receives and processes record data from the `Network` chare array
– `NetData::WriteRecord()` - serializes the record data onto disk

### Halting

After network data and records have been successfully written, control is returned to `Main::FiniSim()` which exits STACS.
