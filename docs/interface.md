# Interface

The interface between the simulation and external devices is handled through YARP.
When STACS is compiled with YARP some additional elements in the programmatic flow become relevant.
By default, there is a method that enables external, interactive control of the simulation by the user.
If applicable, network models that utilize YARP may also be used during simulation to provide external input and output capabilities to a variety of external devices to and from the simulation, respectively.

For STACS to operate using YARP, a YARP server on the network needs to be running.
This may be invoked by calling `yarp server` from the command line.

## STACS with YARP

### Initialization

During initialization, additions to the programmatic flow are as follows:

* In `Main::Main(CkArgMsg *msg)`, during configuration:
  * `Main::BuildVtxDist()` builds an `mVtxDist` message, which provides information about the distribution of vertices across the network partitions
  * An additional chare constructor is called, `StreamRPC::StreamRPC(mVtxDist *msg)`, which enables communication through YARP
* YARP is initialized on each computing unit during the `NetData` constructor call
* Once YARP is initialized, `Main::InitSim()` calls `StreamRPC::Open(CProxy_Network cpnet)` which opens up the RPC port for communication with a reference to the `Network` chare array
* In `Network::LoadNetwork(mPart *msg)`, when the network models are being loaded, any models with YARP ports have them opened for communication

### Simulation

During the simulation, additions to the programmatic flow are as follows:

* RPC messages may be sent to the simulation on the `/stacs/rpc` port (see below)
  * `RPCReader::read(yarp::os::ConnectionReader &connection)` receives any RPC messages through YARP
  * `RPCReader::BuildRPCMsg(idx_t command, yarp::os::Bottle message)` builds an mRPC message, which then can be sent through Charm++
  * `Network::RPCMsg(mRPC *msg)` receives mRPC messages and processes them accordingly on the neural network
* During the main simulation loop, `iter == synciter` is an additional condition that is checked during the control flow
  * If the simulation is not paused or is told to unpause, `StreamRPC::RPCSync()` acts as a coordination point after an `mRPC` message is processed, and control is returned to the main simulation loop
  * If the simulation is paused or is told to pause, `StreamRPC::RPCPause()` acts as a coordination point after an `mRPC` message is processed

### Finalization

During the finalization, additions to the programmatic flow are as follows:

* When the simulation is stopped, `Main::StopSim()` calls `StreamRPC::Close()` which closes the RPC port
* As the neural network is being saved, `Network::SaveNetwork()` closes any YARP ports that have been opened by the network models
* YARP is finalized on each computing unit during the `NetData::SaveNetwork(mPart *msg)` after network data has been serialized

## External Control

With YARP, there is an additional chare in STACS, `StreamRPC`, which is responsible for external, interactive control of the simulation.
Specifically, this sets up an `RPCReader` callback for STACS which accepts messages on the `/stacs/rpc` port.

Even if there are no streams available through network models, STACS provides a way for the simulation to be controlled interactively through the use of remote procedure calls (RPC) through YARP.
This method may be used to send some basic commands for pausing the simulation, stepping through simulation time, and applying manual stimulation commands.
The RPC interface may be invoked by calling `yarp rpc /stacs/rpc` from the command line.

RPC Commands to control simulation are:

* `pause` - pause or unpause the simulation
* `stop` - stop the simulation (final state of the network is written out)
* `check` - checkpoint the simulation (when paused)
* `step <t>` - step the simulation `<t>` _ms_ (when paused)
* `stim <t o a d>` - apply arbitrary stimulation pulses to the network with specified `<t>` targets, `<o>` offsets, `<a>` amplitudes, and `<d>` durations

For the stimulation pulses, the targets may be specified as:

* Individual neurons - both the number of neurons and their indices need to be provided
* Spacial area - both the coordinate of the center point and a radius need to be provided

Each pulse in a stimulation command is parameterized as:

* Offset - time in ms from the time when a stimulation command is given that the pulse starts
* Amplitude - amplitude (unitless) of the pulse that will be received by the targeted network
models
* Duration - time in _ms_ from the time when the pulse starts that the pulse is active

## Streams

YARP streams refer to network models that open ports to either receive or send data from the network simulation online through YARP.
In general, streams that accept more standard YARP signals, such as for audio, will remain better suited to varying devices due to their increased modularity.

A YARP port is opened by the network models that get initialized by the `Network` chare array, which then sets up a callback function, `onRead`, for when data is available on that port.
These ports are decoupled, according to the YARP philosophy, such that the input stream may stop and start at any time without the simulation stalling.
When streaming data from real-world devices, care should be taken so that the complexity of the neural network matches the computing resources in order to run in real-time.
For prerecorded data, all that is needed is for the rate of the streams to be comparable to that of the simulation.

### Audio

The `yarpaudio` model supports streaming audio data through the YARP protocol.
Configuration of this model includes how that audio data will be processed, transforming from a time domain signal into the frequency domain, prior to sending real valued outputs to the simulator.
