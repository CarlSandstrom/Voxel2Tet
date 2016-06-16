#ifndef TIMESTAMP_H
#define TIMESTAMP_H

#include <vector>
#include <string>
#include <ctime>

namespace voxel2tet
{

/**
 * @brief Simple timer to store times between StartTimer and StopTimer members being called
 */
class TimeStamp
{
private:
    std::vector<std::pair<double, std::string>> Table;
    clock_t t0;
    clock_t tend;
    clock_t inittime;
    std::string CurrentMessage;
public:
    TimeStamp();

    /**
     * @brief Clears list and starts first timer;
     * @param Message Message associated with initialization
     */
    void Initialize(std::string Message);

    /**
     * @brief Starts the timer
     * @param Message Message to this item
     */
    void StartTimer(std::string Message);

    /**
     * @brief Stops current timer and stores message and time
     */
    void StopTimer();

    /**
     * @brief Returns time passed between initialization and now
     * @return
     */
    int GiveTotalTime();

    /**
     * @brief Prints table of times elapsed and message for each item
     */
    void PrintTable();

    /**
     * @brief Returns table of times and messages
     * @return Table of times and messages
     */
    std::vector<std::pair<double, std::string>> GetTable();
};

}
#endif // TIMESTAMP_H
