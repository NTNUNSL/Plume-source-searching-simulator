#ifndef SEARCH_ZIGZAG_HPP
#define SEARCH_ZIGZAG_HPP

struct Offset {
    double x;
    double y;
};

// 貪婪飛行法模型
class Zigzag {
private:
    double _direction; // clockwise degree from x axis
    double _drift; // degree from direction
    double _distance; // meters
    char _zigzagStat; // {'r','l'}  current direction

public:
    Zigzag();
    Zigzag(double direction,double drift,double distance);
    
    void trim(double degree);
    void setDirection(double direction){_direction = direction;};
    void setDistance(double distance){_distance = distance;};
    void setDrift(double drift){_drift = drift;};
    Offset genNextOffset();
    char curStat();
    double getDirection(){return _direction;};
    double getDrift(){return _drift;};
    double getDirectionDrifted();
};

#endif  
