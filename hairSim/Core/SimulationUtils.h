#ifndef SIMULATION_UTILS_H
#define SIMULATION_UTILS_H

struct ProxColPointerCompare
{
    bool operator()( const ProximityCollision* lhs, const ProximityCollision* rhs ) const
    {
        return *lhs < *rhs;
    }
};

class SimulationUtils
{
public:
      
      // take all the less important stuff out of StrandImplicitManager/Simulation and put it here
    void updateParameters( const SimulationParameters& params );

    void accumulateProxies( origProxys );

protected:

private:
};

#endif
