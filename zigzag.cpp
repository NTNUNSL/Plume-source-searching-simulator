#include "zigzag.hpp"
#include <cmath>

#define DEG2RAD 0.0174532925

Zigzag::Zigzag(){
    _direction = 0;
    _drift = 30;
    _distance = 10;
    _zigzagStat = 'l';
}

Zigzag::Zigzag(double direction,double drift,double distance){
    _direction = direction;
    _drift = drift;
    _distance = distance;
    _zigzagStat = 'l';
}

void Zigzag::trim(double degree){
    _direction = std::fmod((_direction+degree),360);
}

// 產生鋸齒移動軌跡的下一步
Offset Zigzag::genNextOffset(){
    Offset o;
    double goDir;
    switch(_zigzagStat){
        case 'r':
            _zigzagStat = 'l';
            goDir = _direction - _drift;
            break;
        case 'l':
            _zigzagStat = 'r';
            goDir = _direction + _drift;
            break;
    }
    o.x = _distance * cos((90-goDir) * DEG2RAD);
    o.y = _distance * sin((90-goDir) * DEG2RAD);
    return o;
}

char Zigzag::curStat(){
    return _zigzagStat;
}

// 加上偏移角
double Zigzag::getDirectionDrifted(){
    double s;
    if(_zigzagStat=='r'){
        s = std::fmod(_direction+_drift,360);
    }else{
        s = std::fmod(_direction-_drift,360);
    }
    return s;
}
