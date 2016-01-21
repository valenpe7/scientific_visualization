#ifndef _NODE_
#define _NODE_

#include <iostream>
#include <array>

class node {
	friend class cell;
	friend class weather_data;
public:
	node() = default;
	~node() = default;
private:
	int id;
	std::array<double, 3> values;
};

#endif
