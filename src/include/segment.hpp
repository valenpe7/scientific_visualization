#ifndef _SEGMENT_
#define _SEGMENT_

#include <iostream>
#include <array>

class segment {
	friend class cell;
	friend class weather_data;
public:
	segment() = default;
	~segment() = default;
private:
	std::array<double, 2> start;
	std::array<double, 2> end;
};

#endif
