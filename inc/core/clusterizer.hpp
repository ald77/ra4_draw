#ifndef H_CLUSTERIZER
#define H_CLUSTERIZER

#include <list>
#include <set>
#include <vector>
#include <ostream>
#include <random>

#include "TH2D.h"
#include "TGraph.h"

namespace Clustering{
  class Point{
  public:
    Point() = default;
    Point(float x, float y, float w);
    
    float x_, y_, w_;

    bool operator<(const Point &other) const;
  };

  float WeightedDistance(const Point &a, const Point &b);

  class Node : public Point{
  public:
    Node(float x, float y, float z);
    Node(const Point &p);

    bool operator<(const Node &other) const;

    float dist_to_neighbor_;
    std::list<Node>::iterator neighbor_;
    std::vector<std::list<Node>::iterator> neighbor_of_;
  };

  class Clusterizer{
  public:
    explicit Clusterizer(const TH2D &hist_template,
                         long max_points = -1);

    void AddPoint(float x, float y, float w);
    
    void SetPoints(const std::vector<Point> &points);
    void SetPoints(const TH2D &h);

    TH2D GetHistogram(double luminosity) const;
    TGraph GetGraph(double luminosity, bool keep_in_frame = true) const;

  private:
    long max_points_;
    bool hist_mode_;
    TH2D hist_;
    std::vector<Point> orig_points_;
    mutable std::list<Node> nodes_;
    mutable std::vector<Point> final_points_;
    mutable float clustered_lumi_;

    static std::mt19937_64 prng_;
    static std::uniform_real_distribution<float> urd_;

    void InsertPoint(float x, float y, float w) const;
    void InsertPoint(const Point &p) const;
    void RemovePoint(std::list<Node>::iterator node) const;
    
    std::list<Node>::iterator NearestNeighbors() const;
    
    static void Link(std::list<Node>::iterator node,
                    std::list<Node>::iterator neighbor,
                    float dist);

    void EmptyHistogram();
    void ConvertToHist();

    void Cluster(double luminosity) const;
    void SetupNodes(double luminosity) const;
    void MergeNodes() const;
    void MergeNodes(std::list<Node>::iterator a,
                    std::list<Node>::iterator b) const;
    void SplitNode() const;
  };
}

std::ostream & operator<<(std::ostream &stream, const Clustering::Point &p);
std::ostream & operator<<(std::ostream &stream, const Clustering::Node &n);

#endif
