#include "include/cell.hpp"

segment cell::top_bottom() {
	segment segment;
	segment.start = {
		(this->nodes[2]->values[0] + this->nodes[3]->values[0]) / 2.0,
		(this->nodes[2]->values[1])
	};
	segment.end = {
		(this->nodes[0]->values[0] + this->nodes[1]->values[0]) / 2.0,
		(this->nodes[0]->values[1])
	};
	return segment;
}

segment cell::top_left() {
	segment segment;
	segment.start = {
		(this->nodes[0]->values[0] + this->nodes[1]->values[0]) / 2.0,
		(this->nodes[0]->values[1])
	};
	segment.end = {
		(this->nodes[0]->values[0]),
		(this->nodes[0]->values[1] + this->nodes[3]->values[1]) / 2.0
	};
	return segment;
}

segment cell::top_right() {
	segment segment;
	segment.start = {
		(this->nodes[1]->values[0]),
		(this->nodes[1]->values[1] + this->nodes[2]->values[1]) / 2.0
	};
	segment.end = {
		(this->nodes[0]->values[0] + this->nodes[1]->values[0]) / 2.0,
		(this->nodes[0]->values[1])
	};
	return segment;
}

segment cell::bottom_left() {
	segment segment;
	segment.start = {
		(this->nodes[2]->values[0] + this->nodes[3]->values[0]) / 2.0,
		(this->nodes[2]->values[1])
	};
	segment.end = {
		(this->nodes[0]->values[0]),
		(this->nodes[0]->values[1] + this->nodes[3]->values[1]) / 2.0
	};
	return segment;
}

segment cell::bottom_right() {
	segment segment;
	segment.start = {
		(this->nodes[2]->values[0] + this->nodes[3]->values[0]) / 2.0,
		(this->nodes[2]->values[1])
	};
	segment.end = {
		(this->nodes[1]->values[0]),
		(this->nodes[1]->values[1] + this->nodes[2]->values[1]) / 2.0
	};
	return segment;
}

segment cell::left_right() {
	segment segment;
	segment.start = {
		(this->nodes[0]->values[0]),
		(this->nodes[0]->values[1] + this->nodes[3]->values[1]) / 2.0
	};
	segment.end = {
		(this->nodes[1]->values[0]),
		(this->nodes[1]->values[1] + this->nodes[2]->values[1]) / 2.0
	};
	return segment;
}

cell& cell::compute_average() {
	double sum = 0.0;
	for(std::size_t k = 0; k < 4; ++k) {
		sum += this->nodes[k]->values[2];
	}
	this->average = sum / 4.0;
	return (*this);
}