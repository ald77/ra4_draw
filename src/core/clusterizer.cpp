#include "core/clusterizer.hpp"

#include <cmath>

#include <tuple>
#include <array>
#include <random>

#include "core/utilities.hpp"

using namespace std;
using namespace Clustering;

namespace{
  mt19937_64 InitializePRNG(){
    array<int, 128> sd;
    random_device r;
    generate_n(sd.begin(), sd.size(), ref(r));
    seed_seq ss(begin(sd), end(sd));
    return mt19937_64(ss);
  }
}

Point::Point(float x, float y, float w):
  x_(x),
  y_(y),
  w_(w){  
  }

bool Point::operator<(const Point &other) const{
  return make_tuple(x_, y_, w_)<make_tuple(other.x_, other.y_, other.w_);
}

float Clustering::WeightedDistance(const Point &a, const Point &b){
  float dx = a.x_ - b.x_;
  float dy = a.y_ - b.y_;
  return a.w_*b.w_*(dx*dx+dy*dy)/(a.w_+b.w_);
}

Node::Node(float x, float y, float w):
  Point(x, y, w),
  dist_to_neighbor_(-1.),
  neighbor_(),
  neighbor_of_(){
}

Node::Node(const Point &p):
  Point(p),
  dist_to_neighbor_(-1.),
  neighbor_(),
  neighbor_of_(){
}

bool Node::operator<(const Node &other) const{
  return make_tuple(x_, y_, w_, dist_to_neighbor_)<make_tuple(other.x_, other.y_, other.w_, other.dist_to_neighbor_);
}

mt19937_64 Clusterizer::prng_ = InitializePRNG();
uniform_real_distribution<float> Clusterizer::urd_(0., 1.);

Clusterizer::Clusterizer(const TH2D &hist_template, long max_points):
  max_points_(max_points),
  hist_mode_(max_points == 0),
  hist_(hist_template),
  orig_points_(),
  nodes_(),
  final_points_(),
  clustered_lumi_(-1.){
  if(max_points_ >= 0 && max_points_ < hist_.GetNcells()){
    max_points_ = hist_.GetNcells();
  }
}

void Clusterizer::AddPoint(float x, float y, float w){
  clustered_lumi_ = -1.;
  if(orig_points_.size() >= static_cast<size_t>(max_points_)
     && max_points_ >= 0){
    hist_mode_ = true;
    orig_points_.clear();
  }
  hist_.Fill(x, y, w);
  if(!hist_mode_){
    orig_points_.emplace_back(x, y, w);
  }
}

void Clusterizer::SetPoints(const vector<Point> &points){
  clustered_lumi_ = -1.;
  EmptyHistogram();
  if(points.size() > static_cast<size_t>(max_points_) && max_points_ >= 0){
    hist_mode_ = true;
    orig_points_.clear();
  }else{
    hist_mode_ = false;
    orig_points_ = points;
  }
  for(const auto &p: points){
    hist_.Fill(p.x_, p.y_, p.w_);
  }
}

void Clusterizer::SetPoints(const TH2D &h){
  clustered_lumi_ = -1.;
  EmptyHistogram();
  hist_mode_ = true;
  hist_ = h;
  if(max_points_ >= 0 && max_points_ < hist_.GetNcells()){
    max_points_ = hist_.GetNcells();
  }
}

TH2D Clusterizer::GetHistogram(double luminosity) const{
  TH2D h = hist_;
  h.Scale(luminosity);
  return h;
}

TGraph Clusterizer::GetGraph(double luminosity, bool keep_in_frame) const{
  Cluster(luminosity);
  float xmin = hist_.GetXaxis()->GetBinLowEdge(1);
  float xmax = hist_.GetXaxis()->GetBinUpEdge(hist_.GetNbinsX());
  float dx = 0.0001*(xmax-xmin);
  xmin += dx;
  xmax -= dx;
  float ymin = hist_.GetYaxis()->GetBinLowEdge(1);
  float ymax = hist_.GetYaxis()->GetBinUpEdge(hist_.GetNbinsY());
  float dy = 0.0001*(ymax-ymin);
  ymin += dy;
  ymax -= dy;
  
  TGraph g(final_points_.size());
  g.SetMarkerStyle(hist_.GetMarkerStyle());
  g.SetMarkerColor(hist_.GetMarkerColor());
  g.SetMarkerSize(hist_.GetMarkerSize());
  g.SetTitle(FullTitle(hist_).c_str());

  for(size_t i = 0; i < final_points_.size(); ++i){
    const Point &p = final_points_.at(i);
    float x = p.x_;
    float y = p.y_;
    if(keep_in_frame){
      if(x < xmin) x = xmin;
      if(x > xmax) x = xmax;
      if(y < xmin) y = ymin;
      if(y > ymax) y = ymax;
    }
    g.SetPoint(i, x, y);
  }

  return g;
}

void Clusterizer::InsertPoint(float x, float y, float w) const{
  InsertPoint(Point(x, y, w));
}

void Clusterizer::InsertPoint(const Point &p) const{
  nodes_.emplace_front(p);
  if(nodes_.size() == 2){
    auto first = nodes_.begin();
    auto second = first;
    ++second;
    float dist = WeightedDistance(*first, *second);
    Link(first, second, dist);
    Link(second, first, dist);
  }else if(nodes_.size() == 1){
    return;
  }
  
  list<Node>::iterator this_node = nodes_.begin();
  auto the_end = nodes_.end();
  for(auto it = ++nodes_.begin(); it != the_end; ++it){
    float dist = WeightedDistance(nodes_.front(), *it);
    
    if(dist < this_node->dist_to_neighbor_ || this_node->dist_to_neighbor_ < 0.){
      Link(this_node, it, dist);
    }
    
    if(dist < it->dist_to_neighbor_){
      Link(it, this_node, dist);
    }
  }
}

void Clusterizer::RemovePoint(list<Node>::iterator node) const{
  if(nodes_.size() == 2){
    nodes_.erase(node);
    nodes_.front().dist_to_neighbor_ = -1.;
    nodes_.front().neighbor_ = list<Node>::iterator();
    nodes_.front().neighbor_of_.clear();
    return;
  }else if(nodes_.size() == 1){
    nodes_.clear();
    return;
  }
  
  auto pos = find(node->neighbor_->neighbor_of_.begin(), node->neighbor_->neighbor_of_.end(), node);
  node->neighbor_->neighbor_of_.erase(pos);
  
  auto invalidated = node->neighbor_of_;

  nodes_.erase(node);

  //Set invalid distance on nodes with unknown neighbor
  for(auto &bad_node: invalidated){
    bad_node->dist_to_neighbor_ = -1.;
  }

  //Recompute neighbor of all nodes which had neighbor removed
  for(auto new_neighbor = nodes_.begin(); new_neighbor != nodes_.end(); ++new_neighbor){
    for(auto &bad_node: invalidated){
      if(bad_node == new_neighbor) continue;
      float dist = WeightedDistance(*bad_node, *new_neighbor);
      if(dist < bad_node->dist_to_neighbor_ || bad_node->dist_to_neighbor_ < 0.){
        Link(bad_node, new_neighbor, dist);
      }
    }
  }
}

list<Node>::iterator Clusterizer::NearestNeighbors() const{
  float min_dist = -1., min_dist_high_weight = -1.;
  list<Node>::iterator best_node = nodes_.begin(), best_node_high_weight = nodes_.begin();
  for(auto node = nodes_.begin(); node != nodes_.end(); ++node){
    if(node->dist_to_neighbor_ < min_dist || min_dist < 0.){
      best_node = node;
      min_dist = node->dist_to_neighbor_;
    }
    if((node->dist_to_neighbor_ < min_dist_high_weight || min_dist_high_weight < 0.)
       && node->w_ > 1. && node->neighbor_->w_ < 1.){
      best_node_high_weight = node;
      min_dist_high_weight = node->dist_to_neighbor_;
    }
  }
  if(min_dist < 0.) ERROR("Could not find neighboring points.");
  if(min_dist_high_weight > 0.){
    return best_node_high_weight;
  }else{
    return best_node;
  }
}

void Clusterizer::Link(list<Node>::iterator node,
                       list<Node>::iterator neighbor,
                       float dist){
  //Remove backlink from previous neighbor if necessary
  if(node->dist_to_neighbor_ >= 0.){
    auto pos = find(node->neighbor_->neighbor_of_.begin(), node->neighbor_->neighbor_of_.end(), node);
    node->neighbor_->neighbor_of_.erase(pos);
  }

  //Set up new link
  node->dist_to_neighbor_ = dist;
  node->neighbor_ = neighbor;
  neighbor->neighbor_of_.push_back(node);
}

void Clusterizer::EmptyHistogram(){
  for(int i = 0; i < hist_.GetNcells(); ++i){
    hist_.SetBinContent(i, 0.);
    hist_.SetBinError(i, 0.);
  }
  hist_.SetEntries(0.);
}

void Clusterizer::Cluster(double luminosity) const{
  if(luminosity == clustered_lumi_) return;

  SetupNodes(luminosity);
  MergeNodes();
  
  clustered_lumi_ = luminosity;
}

void Clusterizer::SetupNodes(double luminosity) const{
  nodes_.clear();
  final_points_.clear();
  
  if(hist_mode_){
    int nx = hist_.GetNbinsX();
    int ny = hist_.GetNbinsY();
    float xmin = hist_.GetXaxis()->GetBinLowEdge(1);
    float xmax = hist_.GetXaxis()->GetBinUpEdge(hist_.GetNbinsX());
    float dx = 0.5*(xmax-xmin);
    float ymin = hist_.GetYaxis()->GetBinLowEdge(1);
    float ymax = hist_.GetYaxis()->GetBinUpEdge(hist_.GetNbinsY());
    float dy = 0.5*(ymax-ymin);
    for(int ix = 0; ix <= nx+1; ++ix){
      float xlow = (ix <= 0) ? (xmin-dx)
        : (ix > hist_.GetNbinsX()) ? (xmax+dx)
        : hist_.GetXaxis()->GetBinLowEdge(ix);
      float xhigh = (ix <= 0) ? (xmin-dx)
        : (ix > hist_.GetNbinsX()) ? (xmax+dx)
        : hist_.GetXaxis()->GetBinUpEdge(ix);
      for(int iy = 0; iy <= ny+1; ++iy){
        float ylow = (iy <= 0) ? (ymin-dy)
          : (iy > hist_.GetNbinsY()) ? (ymax+dy)
          : hist_.GetYaxis()->GetBinLowEdge(iy);
        float yhigh = (iy <= 0) ? (ymin-dy)
          : (iy > hist_.GetNbinsY()) ? (ymax+dy)
          : hist_.GetYaxis()->GetBinUpEdge(iy);
        
        float w = luminosity*hist_.GetBinContent(ix, iy);
        while(w > 0.){
          float x = xlow + urd_(prng_)*(xhigh-xlow);
          float y = ylow + urd_(prng_)*(yhigh-ylow);
          if(w >= 1.){
            final_points_.emplace_back(x, y, 1.);
            w -= 1.;
          }else{
            InsertPoint(x, y, w);
            w = 0.;
          }
        }
      }
    }
  }else{
    for(const auto &p: orig_points_){
      float w = luminosity * p.w_;
      if(w <= 0.) continue;
      if(w == 1.){
        final_points_.emplace_back(p.x_, p.y_, w);
      }else{
        InsertPoint(p.x_, p.y_, w);
      }
    }
  }
}

void Clusterizer::MergeNodes() const{
  while(nodes_.size()>0){
    while(nodes_.size()>1){
      list<Node>::iterator root_node = NearestNeighbors();
      list<Node>::iterator neighbor = root_node->neighbor_;
      MergeNodes(root_node, neighbor);
    }

    if(nodes_.size()==1){
      if(nodes_.front().w_ > 1.5){
        SplitNode();
      }else if(nodes_.front().w_ >= 0.5){
        final_points_.push_back(static_cast<Point>(nodes_.front()));
        nodes_.clear();
      }else{
        nodes_.clear();
      }
    }
  }
}

void Clusterizer::MergeNodes(list<Node>::iterator a,
                             list<Node>::iterator b) const{
  if(a->w_ < b->w_){
    //Make sure node "A" has higher weight
    MergeNodes(b, a);
    return;
  }

  if(a->w_ + b->w_ <= 1.){
    //Merge two points into one
    Point c((a->w_*a->x_+b->w_*b->x_)/(a->w_+b->w_),
            (a->w_*a->y_+b->w_*b->y_)/(a->w_+b->w_),
            a->w_+b->w_);
    RemovePoint(a);
    RemovePoint(b);
    if(c.w_ == 1.){
      final_points_.push_back(c);
    }else{
      InsertPoint(c);
    }
  }else{
    //Partition so one point has weight exactly 1
    float sumw = a->w_ + b->w_;
    float summ1 = sumw - 1.;
    float rt = sqrt(a->w_*b->w_*summ1);
    
    if(fabs(1.-a->w_) <= fabs(1.-b->w_)){
      //Transfer weight until A has weight exactly 1
      Point c(((a->w_+rt)*a->x_ + (b->w_-rt)*b->x_)/sumw,
              ((a->w_+rt)*a->y_ + (b->w_-rt)*b->y_)/sumw,
              1.);
      Point d(((a->w_*summ1-rt)*a->x_ + (b->w_*summ1+rt)*b->x_)/(sumw*summ1),
              ((a->w_*summ1-rt)*a->y_ + (b->w_*summ1+rt)*b->y_)/(sumw*summ1),
              summ1);
      RemovePoint(a);
      RemovePoint(b);
      final_points_.push_back(c);
      if(d.w_ == 1.){
        final_points_.push_back(d);
      }else{
        InsertPoint(d);
      }
    }else{
      //Transfer weight until B has weight exactly 1
      Point c(((a->w_*summ1+rt)*a->x_ + (b->w_*summ1-rt)*b->x_)/(sumw*summ1),
              ((a->w_*summ1+rt)*a->y_ + (b->w_*summ1-rt)*b->y_)/(sumw*summ1),
              summ1);
      Point d(((a->w_-rt)*a->x_ + (b->w_+rt)*b->x_)/sumw,
              ((a->w_-rt)*a->y_ + (b->w_+rt)*b->y_)/sumw,
              1.);
      RemovePoint(a);
      RemovePoint(b);
      final_points_.push_back(d);
      if(c.w_ == 1.){
        final_points_.push_back(c);
      }else{
        InsertPoint(c);
      }
    }
  }
}

void Clusterizer::SplitNode() const{
  list<Node>::iterator old = nodes_.begin();
  Point a, b;
  if(final_points_.size() > 0){
    float min_dist = -1.;
    size_t best_index = 0;
    for(size_t index = 0; index < final_points_.size(); ++index){
      float dist = WeightedDistance(*old, final_points_.at(index));
      if(dist < min_dist || min_dist < 0.){
        min_dist = dist;
        best_index = index;
      }
    }
    Point &p = final_points_.at(best_index);
    float dx = old->x_ - p.x_;
    float dy = old->y_ - p.y_;
    float scale = 0.25;
    a = Point(old->x_ + scale*dy, old->y_ - scale*dx, 0.5*old->w_);
    b = Point(old->x_ - scale*dy, old->y_ + scale*dx, 0.5*old->w_);
  }else{
    a = Point(old->x_+1., old->y_+1., 0.5*old->w_);
    b = Point(old->x_-1., old->y_-1., 0.5*old->w_);
  }
  RemovePoint(old);
  if(a.w_ != 1.){
    InsertPoint(a);
  }else{
    final_points_.push_back(a);
  }
  if(b.w_ != 1.){
    InsertPoint(b);
  }else{
    final_points_.push_back(b);
  }
}

ostream & operator<<(ostream &stream, const Clustering::Point &p){
  stream << "Point(x=" << p.x_ << ",y=" << p.y_ << ",w=" << p.w_ << ")";
  return stream;
}

ostream & operator<<(ostream &stream, const Clustering::Node &n){
  stream << "Point(x=" << n.x_ << ",y=" << n.y_ << ",w=" << n.w_ << ",d=" << n.dist_to_neighbor_ << ")";
  return stream;
}
