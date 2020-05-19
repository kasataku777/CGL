#include<iostream>
#include<stdio.h>
#include<math.h>

using namespace std;

#define EPS (1e-10)
#define equals(a,b) (fabs((a)-(b))<EPS)

class Point {
public:
	double x, y;


	Point(double x = 0, double y = 0) :x(x), y(y) {}

	Point operator + (Point p) { return Point(x + p.x, y + p.y); }
	Point operator - (Point p) { return Point(x - p.x, y - p.y); }
	Point operator * (double a) { return Point(a * x, a * y); }
	Point operator / (double a) { return Point(x / a, y / a); }

	double abs() { return sqrt(norm()); }
	double norm() { return x * x + y * y; };

	bool operator < (const Point& p)const {
		return x != p.x ? x < p.x : y < p.y;
	}

	bool operator ==(const Point& p)const {
		return fabs(x - p.x) < EPS && fabs(y - p.y) < EPS;
	}

};

typedef Point Vector;

double dot(Vector a, Vector b) {
	return a.x * b.x + a.y * b.y;
}

double cross(Vector a, Vector b) {
	return a.x * b.y - a.y * b.x;
}

bool isOrthogonal(Vector a, Vector b) {
	return equals(dot(a, b), 0.0);
}

bool isParallel(Vector a, Vector b) {
	return equals(cross(a, b), 0.0);
}
struct Segment {
	Point p1, p2;
	Segment() {}
	Segment(Point p1, Point p2) :p1(p1), p2(p2) {}
};
typedef Segment  Line;

double norm(Vector a) {
	return a.x * a.x + a.y * a.y;
}

Point project(Segment s, Point p) {
	Vector base = s.p2 - s.p1;
	double r = dot(p - s.p1, base) / norm(base);
	return s.p1 + base * r;
}

Point reflect(Segment s, Point p) {
	return p + (project(s, p)-p) * 2.0;
}

int main() {
	double x1, y1, x2, y2, x0, y0;
	cin >> x1 >> y1 >> x2 >> y2;
	int q;
	Line s = Line(Point(x1, y1), Point(x2, y2));

	cin >> q;
	for (int i = 0; i < q; i++) {
		cin >> x0 >> y0;
		Point ans = reflect(s, Point(x0, y0));

		//cout << ans.x << " " << ans.y << endl;
		printf("%.8lf %.8lf\n", ans.x, ans.y);
	}

	return 0;
}