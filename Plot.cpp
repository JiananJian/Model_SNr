#include "Plot.h"
#include "arithmetic.h"

Plot& Plot::xlabel(const char* label) {
	*this << "set xlabel '" << label << "'\n";
	return *this;
}

Plot& Plot::ylabel(const char* label) {
	*this << "set ylabel '" << label << "'\n";
	return *this;
}

Plot& Plot::xrange(double x0, double x1) {
	*this << "set xrange [";
	if (isfinite(x0)) *this << x0;
	*this << ":";
	if (isfinite(x1)) *this << x1;
	*this << "]\n";
	return *this;
}

Plot& Plot::yrange(double x0, double x1) {
	*this << "set yrange [";
	if (isfinite(x0)) *this << x0;
	*this << ":";
	if (isfinite(x1)) *this << x1;
	*this << "]\n";
	return *this;
}

Plot& Plot::key(bool on) {
	if (on)
		*this << "set key\n";
	else
		*this << "unset key\n";
	return *this;
}

Plot& Plot::save(const char* file) {
	*this << "set terminal png\n";
	*this << "set output '" << file << "'\n";
	this->flush(); // *this << "replot\n";
	*this << "unset output\n";
	*this << "unset terminal\n";
	return *this;
}

Plot& Plot::plot(std::vector<double> x0, std::vector<double> y0, const char* title, const char* style) {
	x.push_back(x0);
	y.push_back(y0);
	std::string _title = title;
	std::string _style = style;
	titles.push_back(_title);
	styles.push_back(_style);
	return *this;
}

Plot& Plot::plot(std::vector<double> x, double* y, const char* legend, const char* style) {
	int n = x.size();
	std::vector<double> z(y, y + n);
	return plot(x, z, legend, style);
}

Plot& Plot::plot(std::vector<double> x, double y, const char* legend, const char* style) {
	return plot(x, 0.0 * x + y, legend, style);
}

Plot& Plot::flush() {
	*this << "plot ";
	for (int i = 0; i < x.size(); i++) {
		*this << "'-' with "<< styles[i] <<" title '" << titles[i] << "', ";
	}
	*this << "\n";
	for (int i = 0; i < x.size(); i++) {
		//int j = 1;
		//while (j + 1000 < x.size()) {
		//	std::vector<double> x0(x[i].begin() + j - 1, x[i].begin() + j + 1000+1);
		//	std::vector<double> y0(y[i].begin() + j - 1, y[i].begin() + j + 1000+1);
		//	send1d(boost::make_tuple(x0, y0));
		//	j += 1000;
		//	Gnuplot::flush();
		//}
		//std::vector<double> x0(x[i].begin() + j - 1, x[i].end());
		//std::vector<double> y0(y[i].begin() + j - 1, y[i].end());
		//send1d(boost::make_tuple(x0, y0));
		//Gnuplot::flush();

		send1d(boost::make_tuple(sampling(x[i]), sampling(y[i]))); // limit resolution to accelerate the graphics, but spike peak amplitude may be irregularly reduced
	}
	Gnuplot::flush();
	return *this;
}

Plot& Plot::clear() {
	x.clear();
	y.clear();
	titles.clear();
	styles.clear();
	*this << "clear\n";
	*this << "reset\n";
	return *this;
}

//Plot& Plot::operator<<(const char* string) {
//	gp << string;
//	return *this;
//}
//
//Plot& Plot::operator<<(std::ostream& (*manipulator)(std::ostream&)) {
//	gp << manipulator;
//	return *this;
//}

// Example : how to use gnuplot 
//void f() {
//	double x[5] = { 1.0, 2.0, 3.0, 4.0, 5.0 };
//	double y[5] = { 1.0, 4.0, 9.0, 16.0, 25.0 };
//
//	Gnuplot gp;
//	gp << "plot '-' with lines title 'x vs y'\n";
//	gp.send1d(boost::make_tuple(x, y));
//}