#include "network.h"

// Using yaml-cpp
#include "yaml-cpp/yaml.h"

extern /*readonly*/ std::string netwkdir;

class StimTimeInput : public ModelTmpl < 84, StimTimeInput > {
    public:
        StimTimeInput() {
            // This model does not require keeping track of state
            // Resize all state lists to zero
            statelist.resize(0);
            sticklist.resize(0);
            auxstate.resize(0);
            auxstick.resize(0);

            // ports
            portlist.resize(1);
            portlist[0] = "input";
         }

        // Simulation functions
        tick_t Step(tick_t tdrift, tick_t diff, std::vector<real_t>& state, std::vector<tick_t>& stick, std::vector<event_t>& events);
        void Jump(const event_t& event, std::vector<std::vector<real_t>>& state, std::vector<std::vector<tick_t>>& stick, const std::vector<auxidx_t>& auxidx) {/*nothing to do*/ }

        // Protocol
        void OpenPorts();
        void ClosePorts();

    private:
        YAML::Node input; // yaml containing externally-supplied spikes
        std::list<idx_t> spike_indices;  // List of neuron indices indicating target of each external spike
        std::list<real_t> spike_times;  // List of time stamps for each external spike
        std::list<tick_t> spike_ticks; // spike times expressed in ticks
        idx_t spidx; // the next spike index to emit
        tick_t sptime; // the next spike time to emit
};

// Simulation step
tick_t StimTimeInput::Step(tick_t tdrift, tick_t tdiff, std::vector<real_t>& state, std::vector<tick_t>& stick, std::vector<event_t>& events) {

    if (spike_indices.empty()) {
        return tdiff;
    }
    // pop off spikes from the spike lists that correspond to this tick
    while (sptime <= tdrift) {
        // create spike event
        event_t event;
        event.diffuse = sptime;
        event.type = EVENT_STIM;
        event.source = REMOTE_EDGE; // Each input spike is sent through a specific edge to a specific neuron; use edge singlecast

        event.index = spidx; // Target of this spike
        event.data = 20;
        events.push_back(event); // Add event to event list
        event.data = -20;
        event.diffuse = sptime + TICKS_PER_MS;
        events.push_back(event);

        //CkPrintf("Emitting external spike: %lu, %lu\n", spidx, sptime);
        // remove spike
        spike_indices.pop_front();
        spike_times.pop_front();
        spike_ticks.pop_front();
        if (!spike_indices.empty()) {
            spidx = spike_indices.front();
            sptime = spike_ticks.front();
        } else {
            return tdiff;
        }
    }

    return tdiff;
}

// Read external spike times
void StimTimeInput::OpenPorts() {
    char ymlfile[100];
    sprintf(ymlfile, "%s/%s", netwkdir.c_str(), portname[0].c_str());
    CkPrintf("Reading external spike times from %s\n", ymlfile);
    try {
        input = YAML::LoadFile(ymlfile);
    } catch (YAML::BadFile& error) {
        CkPrintf("  %s\n", error.what());
    }
    
    try {
        spike_indices = input["index"].as<std::list<idx_t>>();
    } catch (YAML::RepresentationException& error) {
        CkPrintf("  warning: external spike indices not defined\n");
    }
    try {
        spike_times = input["times"].as<std::list<real_t>>();
    } catch (YAML::RepresentationException& error) {
        CkPrintf("  warning: external spike times not defined\n");
    }

    // if number of indices does not equal number of times, clear all lists and emit warning
    if (spike_indices.size() != spike_times.size()) {
        CkPrintf(" warning: indices/times length mismatch\n");
        spike_indices.clear();
        spike_times.clear();
    }
    // build spike ticks by multiplying by TICKS_PER_MS
    //
    std::list<real_t>::iterator it;
    for (it = spike_times.begin(); it != spike_times.end(); ++it) {
        spike_ticks.push_back((tick_t)(*it*TICKS_PER_MS));
    }

    spidx = spike_indices.front();
    sptime = spike_ticks.front();
    CkPrintf(" loaded %lu external spikes\n", spike_indices.size());
}

void StimTimeInput::ClosePorts() {
    // nothing to do here
}
