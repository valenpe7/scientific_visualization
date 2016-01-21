#include "include/weather_data.hpp"
#include "include/lodepng/lodepng.h"
#include "include/eigen/LU"

weather_data::weather_data() {
	std::cout << "creating datastructure for weather data processing" << std::endl;
}

void weather_data::load_data(std::string filename) {
	std::ifstream input_file;
	std::string line;
	std::cout << "loading scattered data: " << filename << std::endl;
	this->name = filename.substr(0, filename.find_last_of("."));
	input_file.open(filename);
	if(input_file.fail()) {
		std::cerr << "error: cannot read input file" << std::endl;
		return;
	}
	this->scattered.size = 0;
	while(std::getline(input_file, line)) {
		++this->scattered.size;
	}
	input_file.clear();
	input_file.seekg(0, std::ios::beg);
	this->scattered.values.resize(this->scattered.size);
	if(input_file.is_open()) {
		for(std::size_t i = 0; i < this->scattered.size; ++i) {
			for(std::size_t j = 0; j < 3; ++j) {
				input_file >> this->scattered.values[i][j];
				input_file.get();
			}
		}
	}
	double min_x = this->scattered.values[0][0];
	double max_x = this->scattered.values[0][0];
	double min_y = this->scattered.values[0][1];
	double max_y = this->scattered.values[0][1];
	for(std::size_t i = 0; i < this->scattered.size; ++i) {
		min_x = std::min(min_x, this->scattered.values[i][0]);
		max_x = std::max(max_x, this->scattered.values[i][0]);
		min_y = std::min(min_y, this->scattered.values[i][1]);
		max_y = std::max(max_y, this->scattered.values[i][1]);
	}
	this->scattered.boundary = {min_x, max_x, min_y, max_y};
	std::cout << "scattered data loaded successfully" << std::endl;
	return;
}

void weather_data::create_grid(std::array<unsigned int, 2> size) {
	if(!size[0] || !size[1]) {
		std::cerr << "error: cannot create grid" << std::endl;
		return;
	}
	std::cout << "creating uniform grid of size: " << size[0] << " x " << size[1] << std::endl;
	this->uniform.size = {
		size[0],
		size[1]
	};
	this->uniform.boundary = this->scattered.boundary;
	this->uniform.spacing = {
		(this->uniform.boundary[1] - this->uniform.boundary[0]) / this->uniform.size[0],
		(this->uniform.boundary[3] - this->uniform.boundary[2]) / this->uniform.size[1]
	};
	this->uniform.nodes.resize(this->uniform.size[1]);
	for(std::size_t i = 0; i < this->uniform.size[1]; ++i) {
		this->uniform.nodes[i].resize(this->uniform.size[0]);
	}
	for(std::size_t i = 0; i < this->uniform.size[1]; ++i) {
		for(std::size_t j = 0; j < this->uniform.size[0]; ++j) {
			this->uniform.nodes[i][j].id = this->uniform.size[0] * i + j;
			this->uniform.nodes[i][j].values[0] = this->uniform.boundary[0] + j * this->uniform.spacing[0];
			this->uniform.nodes[i][j].values[1] = this->uniform.boundary[2] + i * this->uniform.spacing[1];
			this->uniform.nodes[i][j].values[2] = 0.0;
		}
	}
	this->uniform.cells.resize(this->uniform.size[1] - 1);
	for(std::size_t i = 0; i < this->uniform.size[1] - 1; ++i) {
		this->uniform.cells[i].resize(this->uniform.size[0] - 1);
	}
	for(std::size_t i = 0; i < this->uniform.size[1] - 1; ++i) {
		for(std::size_t j = 0; j < this->uniform.size[0] - 1; ++j) {
			this->uniform.cells[i][j].id = (this->uniform.size[0] - 1) * i + j;
			this->uniform.cells[i][j].average = 0.0;
			this->uniform.cells[i][j].RGBA = {0, 0, 0, 255};
			this->uniform.cells[i][j].nodes = {
				&this->uniform.nodes[i][j],
				&this->uniform.nodes[i][j + 1],
				&this->uniform.nodes[i + 1][j + 1],
				&this->uniform.nodes[i + 1][j],
			};
		}
	}
	std::cout << "uniform grid created successfully" << std::endl;
	return;
}

weather_data& weather_data::shepard_interpolation(double p) {
	std::cout << "performing interpolation using the sheppard method with parameter p = " << p << std::endl;
	std::vector<std::vector<double>> w;
	w.resize(this->uniform.size[1]);
	for(std::size_t i = 0; i < this->uniform.size[1]; ++i) {
		w[i].resize(this->uniform.size[0]);
	}
	for(std::size_t i = 0; i < this->uniform.size[1]; ++i) {
		for(std::size_t j = 0; j < this->uniform.size[0]; ++j) {
			w[i][j] = 0.0;
			this->uniform.nodes[i][j].values[2] = 0.0;
			for(std::size_t k = 0; k < this->scattered.size; ++k) {
				w[i][j] += 1.0 / pow(sqrt(pow(this->uniform.nodes[i][j].values[0] - this->scattered.values[k][0], 2)
					+ pow(this->uniform.nodes[i][j].values[1] - this->scattered.values[k][1], 2)), p);
			}
			for(std::size_t k = 0; k < this->scattered.size; ++k) {
				this->uniform.nodes[i][j].values[2] += (1.0 / pow(sqrt(pow(this->uniform.nodes[i][j].values[0]
					- this->scattered.values[k][0], 2) + pow(this->uniform.nodes[i][j].values[1]
					- this->scattered.values[k][1], 2)), p)) * (1.0 / w[i][j]) * this->scattered.values[k][2];
			}
		}
	}
	std::cout << "scattered data successfully interpolated to uniform grid" << std::endl;
	return (*this);
}

weather_data& weather_data::hardy_interpolation(double R) {
	std::cout << "performing interpolation using the hardy multiquadrics with parameter R = " << R << std::endl;
	Eigen::MatrixXd b(this->scattered.size, this->scattered.size);
	Eigen::VectorXd alpha(this->scattered.size);
	Eigen::VectorXd f(this->scattered.size);
	for(std::size_t k = 0; k < this->scattered.size; ++k) {
		f(k) = this->scattered.values[k][2];
		for(std::size_t l = 0; l < this->scattered.size; ++l) {
			b(k, l) = sqrt(R + pow(this->scattered.values[k][0] - this->scattered.values[l][0], 2)
				+ pow(this->scattered.values[k][1] - this->scattered.values[l][1], 2));
		}
	}
	alpha = b.lu().solve(f);
	std::vector<double> h;
	h.resize(this->scattered.size);
	for(std::size_t i = 0; i < this->uniform.size[1]; ++i) {
		for(std::size_t j = 0; j < this->uniform.size[0]; ++j) {
			this->uniform.nodes[i][j].values[2] = 0.0;
			for(std::size_t k = 0; k < this->scattered.size; ++k) {
				h[k] = sqrt(R + pow(this->uniform.nodes[i][j].values[0] - this->scattered.values[k][0], 2)
					+ pow(this->uniform.nodes[i][j].values[1] - this->scattered.values[k][1], 2));
				this->uniform.nodes[i][j].values[2] += alpha(k) * h[k];
			}
		}
	}
	std::cout << "scattered data successfully interpolated to uniform grid" << std::endl;
	return (*this);
}

weather_data& weather_data::compute_colormap(std::vector<std::array<int, 4>> color_table) {
	std::cout << "computing colormap.. " << std::endl;
	double max = 0.0;
	double min = std::numeric_limits<double>::infinity();
	for(std::size_t i = 0; i < this->uniform.size[1] - 1; ++i) {
		for(std::size_t j = 0; j < this->uniform.size[0] - 1; ++j) {
			this->uniform.cells[i][j].compute_average();
			max = std::max(max, this->uniform.cells[i][j].average);
			min = std::min(min, this->uniform.cells[i][j].average);
		}
	}
	unsigned int n = color_table.size();
	double span = (max - min) / n;
	std::cout << "number of colors: " << n << std::endl;
	this->colormap.resize((this->uniform.size[0] - 1) * (this->uniform.size[1] - 1) * 4);
	for(std::size_t i = 0; i < this->uniform.size[1] - 1; ++i) {
		for(std::size_t j = 0; j < this->uniform.size[0] - 1; ++j) {
			for(std::size_t k = 0; k < n; ++k) {
				if(this->uniform.cells[i][j].average <= min + (k + 1) * span && this->uniform.cells[i][j].average >= min + k * span) {
					this->uniform.cells[i][j].RGBA = color_table[k];
				}
			}
			for(std::size_t k = 0; k < 4; ++k) {
				this->colormap[4 * ((this->uniform.size[1] - i - 2) * (this->uniform.size[0] - 1) + j) + k] = this->uniform.cells[i][j].RGBA[k];
			}
		}
	}
	std::cout << "colormap created successfully" << std::endl;
	return (*this);
}

weather_data& weather_data::compute_isolines(unsigned int n) {
	std::cout << "computing isolines.. " << std::endl;
	this->isolines.resize(n);
	std::cout << "number of isolines: " << n << std::endl;
	std::vector<double> thresholds;
	double max = 0.0;
	double min = std::numeric_limits<double>::infinity();
	for(std::size_t i = 0; i < this->uniform.size[1]; ++i) {
		for(std::size_t j = 0; j < this->uniform.size[0]; ++j) {
			max = std::max(max, this->uniform.nodes[i][j].values[2]);
			min = std::min(min, this->uniform.nodes[i][j].values[2]);
		}
	}
	double span = (max - min) / (n + 1);
	for(std::size_t i = 0; i < n; i++) {
		thresholds.push_back(min + i * span);
	}
	unsigned int bits;
	for(size_t l = 0; l < this->isolines.size(); ++l) {
		this->isolines[l].value = thresholds[l];
		for(std::size_t i = 0; i < this->uniform.size[1] - 1; ++i) {
			for(std::size_t j = 0; j < this->uniform.size[0] - 1; ++j) {
				bits = 0;
				for(std::size_t k = 0; k < 4; ++k) {
					if(this->uniform.cells[i][j].nodes[k]->values[2] > this->isolines[l].value) {
						bits |= 1 << (3 - k);
					}
				}
				switch(bits) {
					case 0: case 15:
						break;
					case 1: case 14:
						this->isolines[l].segments.push_back(this->uniform.cells[i][j].bottom_left());
						break;
					case 2: case 13:
						this->isolines[l].segments.push_back(this->uniform.cells[i][j].bottom_right());
						break;
					case 3: case 12:
						this->isolines[l].segments.push_back(this->uniform.cells[i][j].left_right());
						break;
					case 4: case 11:
						this->isolines[l].segments.push_back(this->uniform.cells[i][j].top_right());
						break;
					case 6: case 9:
						this->isolines[l].segments.push_back(this->uniform.cells[i][j].top_bottom());
						break;
					case 7: case 8:
						this->isolines[l].segments.push_back(this->uniform.cells[i][j].top_left());
						break;
					case 5:
						if(this->uniform.cells[i][j].average > this->isolines[l].value) {
							this->isolines[l].segments.push_back(this->uniform.cells[i][j].top_left());
							this->isolines[l].segments.push_back(this->uniform.cells[i][j].bottom_right());
						} else {
							this->isolines[l].segments.push_back(this->uniform.cells[i][j].top_right());
							this->isolines[l].segments.push_back(this->uniform.cells[i][j].bottom_left());
						}
						break;
					case 10:
						if(this->uniform.cells[i][j].average > this->isolines[l].value) {
							this->isolines[l].segments.push_back(this->uniform.cells[i][j].top_right());
							this->isolines[l].segments.push_back(this->uniform.cells[i][j].bottom_left());
						}
						else {
							this->isolines[l].segments.push_back(this->uniform.cells[i][j].top_left());
							this->isolines[l].segments.push_back(this->uniform.cells[i][j].bottom_right());
						}
						break;
					default:
						std::cerr << "error: unknown error" << std::endl;
						break;
				}
			}
		}
	}
	std::cout << "isolines computed successfully" << std::endl;
	return (*this);
}

void weather_data::export_colormap() {
	std::cout << "exporting: " << this->name + "_colormap.kml" << std::endl;
	unsigned error = lodepng::encode(this->name + "_colormap.png", this->colormap, this->uniform.size[0] - 1, this->uniform.size[1] - 1);
	if(error) {
		std::cerr << "encoder error " << error << ": " << lodepng_error_text(error) << std::endl;
		return;
	}
	std::ofstream file(this->name + "_colormap.kml");
	if(file) {
		file
			<< "<?xml version= \"1.0\" encoding= \"UTF-8\"?> " << std::endl
			<< "<kml xmlns=\"http://www.opengis.net/kml/2.2\">" << std::endl
			<< "<GroundOverlay>" << std::endl
			<< "<color>" << "66ffffff" << "</color>" << std::endl
			<< "<name>" << this->name + "_colormap" << "</name>" << std::endl
			<< "<Icon>" << std::endl
			<< "<href>" << this->name + "_colormap.png" << "</href>" << std::endl
			<< "</Icon>" << std::endl
			<< "<LatLonBox>" << std::endl
			<< "<east>" << this->uniform.boundary[0] << "</east>" << std::endl
			<< "<west>" << this->uniform.boundary[1] << "</west>" << std::endl
			<< "<south>" << this->uniform.boundary[2] << "</south>" << std::endl
			<< "<north>" << this->uniform.boundary[3] << "</north>" << std::endl
			<< "<rotation>" << 0.0 << "</rotation>" << std::endl
			<< "</LatLonBox>" << std::endl
			<< "</GroundOverlay>" << std::endl
			<< "</kml>" << std::endl;
		file.close();
	} else {
		std::cerr << "error: unable to save kml file" << std::endl;
		return;
	}
	return;
}

void weather_data::export_isolines() {
	std::cout << "exporting: " << this->name + "_isolines.kml" << std::endl;
	std::ofstream file(this->name + "_isolines.kml");
	if(file) {
		file
			<< "<?xml version= \"1.0\" encoding= \"UTF-8\"?> " << std::endl
			<< "<kml xmlns=\"http://www.opengis.net/kml/2.2\">" << std::endl
			<< "<Document>" << std::endl
			<< "<name>" << this->name + "_isolines" << "</name>" << std::endl
			<< "<Style id=\"isolines\">" << std::endl
			<< "<LineStyle>" << std::endl
			<< "<color>" << "66ffffff" << "</color>" << std::endl
			<< "<width>" << 2 << "</width>" << std::endl
			<< "</LineStyle>" << std::endl
			<< "</Style>" << std::endl;
		for(std::size_t i = 0; i < this->isolines.size(); ++i) {
			for(std::size_t j = 0; j < this->isolines[i].segments.size(); ++j) {
				file
					<< std::endl
					<< "<Placemark>" << std::endl
					<< "<styleUrl>" << "#isolines" << "</styleUrl>" << std::endl
					<< "<LineString>" << std::endl
					<< "<coordinates>" << std::endl
					<< this->isolines[i].segments[j].start[0] << ","
					<< this->isolines[i].segments[j].start[1]
					<< std::endl
					<< this->isolines[i].segments[j].end[0] << ","
					<< this->isolines[i].segments[j].end[1]
					<< std::endl
					<< "</coordinates>" << std::endl
					<< "</LineString>" << std::endl
					<< "</Placemark>" << std::endl
					<< std::endl;
			}
		}
		file
			<< "</Document>" << std::endl
			<< "</kml>" << std::endl;
		file.close();
	} else {
		std::cerr << "error: unable to save kml file" << std::endl;
		return;
	}
	return;
}

weather_data::~weather_data() {
	std::cout << "deleting datastructure for weather data processing" << std::endl;
}