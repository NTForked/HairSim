
#ifndef CONFIG_H
#define CONFIG_H


std::cout << "config should not exist, these should be options/parameters in the scenes, sim params files, that get printed out" << endl;

// PARAMETERS
//---------------------------------- CD Loop
const bool loopCD = false;
//---------------------------------- Velocity Condition
const bool velCondition = false;
//---------------------------------- Thickness
const bool thickness = false;
//---------------------------------- Penetration Response Vec to influence RelVel
const bool penetrationResponse = true;
//---------------------------------- Threshold below which redundant collisions are discarded
const double abscissaReplacementTolerance = 0.0; //1.4; // 1.4 ?
//---------------------------------- Ignore Self-Collisions
const bool rejectSelfCollisions = true;

static bool nonLinearCallbackBogus = false;
static bool hLoop = false;
static bool trackGeometricRelations = true;
static bool penaltyOnce = true;
static bool penaltyAfter = true;

#endif
