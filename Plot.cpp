#include "Plot.h"

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
	titles.push_back(title);
	styles.push_back(style);
	return *this;
}

Plot& Plot::plot(std::vector<double> x0, double y0, const char* title, const char* style) {
	x.push_back(x0);
	y.push_back(0.0 * x0 + y0);
	titles.push_back(title);
	styles.push_back(style);
	return *this;
}

Plot& Plot::flush() {
	*this << "plot ";
	for (int i = 0; i < x.size(); i++) {
		*this << "'-' with "<< styles[i] <<" title '" << titles[i] << "', ";
	}
	*this << "\n";
	for (int i = 0; i < x.size(); i++) {
		send1d(boost::make_tuple(x[i], y[i]));
	}
	Gnuplot::flush();
	return *this;
}

Plot& Plot::clear() {
	x.clear();
	y.clear();
	titles.clear();
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