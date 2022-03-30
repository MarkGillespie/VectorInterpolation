#include "utils.h"

double cot(double x) { return cos(x) / sin(x); }

verbose_runtime_error::verbose_runtime_error(const std::string& arg,
                                             const char* file, int line)
    : std::runtime_error(arg) {

    std::string filePath(file);

    // stolen from polyscope/utilities.cpp
    size_t startInd = 0;
    for (std::string sep : {"/", "\\"}) {
        size_t pos = filePath.rfind(sep);
        if (pos != std::string::npos) {
            startInd = std::max(startInd, pos + 1);
        }
    }

    std::string niceName = filePath.substr(startInd, filePath.size());

    std::ostringstream o;
    o << arg << " At " << niceName << ":" << line;
    msg = o.str();
}

Eigen::Vector3d toEigen(const Vector3& v) {
    Eigen::Vector3d ret;
    ret << v.x, v.y, v.z;
    return ret;
}

Vector3 fromEigen(const Eigen::Vector3d& v) {
    return Vector3{v[0], v[1], v[2]};
}

// Returns vw^T
Eigen::Matrix3d outer(const Vector3& v, const Vector3& w) {
    return toEigen(v) * toEigen(w).transpose();
}


string getFilename(string filePath, bool withExtension, char seperator) {
    // Get last dot position
    size_t dotPos = filePath.rfind('.');
    size_t sepPos = filePath.rfind(seperator);

    if (sepPos != string::npos) {
        return filePath.substr(
            sepPos + 1,
            filePath.size() -
                (withExtension || dotPos != std::string::npos ? 1 : dotPos));
    }
    return "";
}
