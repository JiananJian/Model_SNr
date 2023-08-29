#pragma once

#include "gnuplot-iostream.h"
#include "io_data.h"

/**
 * @brief 
 */
class Plot : public Gnuplot {
public:
	//Gnuplot gp;
	//Plot& operator<<(const char* string);
	//Plot& operator<<(std::ostream& (*manipulator)(std::ostream&)); // e.g. std::endl;
	std::vector<std::vector<double>> x, y;
	std::vector<const char*> titles, styles;
	Plot& xlabel(const char* label);
	Plot& ylabel(const char* label);
	Plot& xrange(double x0, double x1);
	Plot& yrange(double x0, double x1);
	Plot& key(bool on);
	Plot& save(const char* file);
	Plot& plot(std::vector<double> x, std::vector<double> y, const char* legend, const char* style = "lines");
	Plot& plot(std::vector<double> x, double y, const char* legend, const char* style = "lines");
	Plot& flush();
	Plot& clear();
};
