#include<iostream>
#include<stdio.h>
#include<math.h>
#include<algorithm>
#include<string>
#include<map>

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

class Circle {
public:
	Point c;
	double r;
	Circle(Point c=Point(),double r=0.0):c(c),r(r){}
};

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
double abs(Vector a) {
	return sqrt(norm(a));
}

Point project(Segment s, Point p) {
	Vector base = s.p2 - s.p1;
	double r = dot(p - s.p1, base) / norm(base);
	return s.p1 + base * r;
}

Point reflect(Segment s, Point p) {
	return p + (project(s, p) - p) * 2.0;
}
double getDistance(Point a, Point b) {
	return abs(a - b);
}

double getDistanceLP(Line l, Point p) {
	return fabs(cross(l.p2 - l.p1, p - l.p1) / abs(l.p2 - l.p1));
}

double getDistanceSP(Segment s, Point p) {
	if (dot(s.p2 - s.p1, p - s.p1) < 0.0)return abs(p - s.p1);
	if (dot(s.p1 - s.p2, p - s.p2) < 0.0)return abs(p - s.p2);
	return getDistanceLP(s, p);

}


static const int COUNTER_CLOCKWISE = 1;
static const int CLOCKWISE = -1;
static const int ONLINE_BACK = 2;
static const int ONLINE_FRONT = -2;
static const int ON_SEGMENT = 0;

int ccw(Point p0, Point p1, Point p2) {
	Vector a = p1 - p0;
	Vector b = p2 - p0;
	if (cross(a, b) > EPS)return COUNTER_CLOCKWISE;
	if (cross(a, b) < -EPS) return CLOCKWISE;
	if (dot(a, b) < -EPS)return ONLINE_BACK;
	if (a.norm() < b.norm()) return ONLINE_FRONT;

	return ON_SEGMENT;


}
bool intersect(Point p1, Point p2, Point p3, Point p4) {
	return (ccw(p1, p2, p3) * ccw(p1, p2, p4) <= 0 && ccw(p3, p4, p1) * ccw(p3, p4, p2) <= 0);
}
bool intersect(Segment s1, Segment s2) {
	return intersect(s1.p1, s1.p2, s2.p1, s2.p2);
}

double getDistance(Segment s1, Segment s2) {
	if (intersect(s1, s2))return 0.0;
	return min(min(getDistanceSP(s1, s2.p1), getDistanceSP(s1, s2.p2)),
		min(getDistanceSP(s2, s1.p1), getDistanceSP(s2, s1.p2)));
}

Point getCrossPoint(Segment s1, Segment s2) {
	Vector base = s2.p2 - s2.p1;
	double d1 = fabs(cross(base, s1.p1 - s2.p1));
	double d2 = fabs(cross(base, s1.p2 - s2.p1));
	double t = d1 / (d1 + d2);
	return s1.p1 + (s1.p2 - s1.p1) * t;
}


pair<Point, Point> getCrossPoints(Circle c, Line l) {
	//assert(intersect(c, l));
	Vector pr = project(l,c.c);
	Vector e = (l.p2 - l.p1) / abs(l.p2 - l.p1);
	double base = sqrt(c.r * c.r - norm(pr - c.c));
	return make_pair(pr + e * base, pr - e * base);
}

double arg(Vector p) { return atan2(p.y, p.x); }
Vector polar(double a, double r) { return Point(cos(r) * a, sin(r) * a); }

pair<Point, Point>getCrossPoints(Circle c1, Circle c2) {
	//assert(intersect(c1,c2);
	double d = abs(c1.c - c2.c);
	double a = acos((c1.r * c1.r + d * d - c2.r * c2.r) / (2 * c1.r * d));
	double t = arg(c2.c - c1.c);
	return make_pair(c1.c + polar(c1.r, t + a), c1.c + polar(c1.r, t - a));



}

//int main() {
//	double x, y, r,x1,y1,x2,y2;
//	int q;
//	pair<Point, Point> p;
//
//	cin >> x >> y >> r;
//	cin >> q;
//
//	for (int i = 0; i < q; i++) {
//		cin >> x1 >> y1 >> x2 >> y2;
//		p= getCrossPoints(Circle(Point(x, y), r), Line(Point(x1, y1), Point(x2, y2)));
//		if (p.second < p.first)swap(p.first, p.second);
//
//		printf("%.8lf %.8lf %.8lf %.8lf\n", p.first.x, p.first.y, p.second.x, p.second.y);
//	}
//
//
//	return 0;
//}

int main() {
	int x1, x2, r1, y1, y2, r2;

	cin >> x1 >> y1 >> r1;
	cin >> x2 >> y2 >> r2;

	pair<Point, Point>p;
	p=getCrossPoints(Circle(Point(x1, y1), r1), Circle(Point(x2, y2), r2));
	if (p.second < p.first)swap(p.first, p.second);
	printf("%.8lf %.8lf %.8lf %.8lf\n", p.first.x, p.first.y, p.second.x, p.second.y);
}