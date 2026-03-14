#ifndef VOXEL2TETEXCEPTION_H
#define VOXEL2TETEXCEPTION_H

#include <stdexcept>
#include <string>

namespace voxel2tet
{
class Voxel2TetException : public std::runtime_error
{
public:
    explicit Voxel2TetException(const std::string &message) : std::runtime_error(message)
    {}
};
}

#endif // VOXEL2TETEXCEPTION_H
