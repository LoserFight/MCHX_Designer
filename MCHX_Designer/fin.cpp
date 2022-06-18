#include "fin.h"

fin::fin(double xFPI, double xfinTickness):FPI(xFPI),finTickness(xfinTickness)
{
	finSpacing = 0.0254 / FPI - finTickness;
}

void fin::geoCalculate()
{
	finSpacing = 0.0254 / FPI - finTickness;
}
