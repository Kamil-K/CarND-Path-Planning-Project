#include <fstream>
#include <math.h>
#include <uWS/uWS.h>
#include <chrono>
#include <iostream>
#include <thread>
#include <vector>
#include "Eigen-3.3/Eigen/Core"
#include "Eigen-3.3/Eigen/QR"
#include "json.hpp"
#include "spline.h"

using namespace std;
using namespace tk;

// for convenience
using json = nlohmann::json;

// For converting back and forth between radians and degrees.
constexpr double pi() { return M_PI; }
double deg2rad(double x) { return x * pi() / 180; }
double rad2deg(double x) { return x * 180 / pi(); }

// Checks if the SocketIO event has JSON data.
// If there is data the JSON object in string format will be returned,
// else the empty string "" will be returned.
string hasData(string s) {
  auto found_null = s.find("null");
  auto b1 = s.find_first_of("[");
  auto b2 = s.find_first_of("}");
  if (found_null != string::npos) {
    return "";
  } else if (b1 != string::npos && b2 != string::npos) {
    return s.substr(b1, b2 - b1 + 2);
  }
  return "";
}

double distance(double x1, double y1, double x2, double y2)
{
	return sqrt((x2-x1)*(x2-x1)+(y2-y1)*(y2-y1));
}
int ClosestWaypoint(double x, double y, const vector<double> &maps_x, const vector<double> &maps_y)
{

	double closestLen = 100000; //large number
	int closestWaypoint = 0;

	for(int i = 0; i < maps_x.size(); i++)
	{
		double map_x = maps_x[i];
		double map_y = maps_y[i];
		double dist = distance(x,y,map_x,map_y);
		if(dist < closestLen)
		{
			closestLen = dist;
			closestWaypoint = i;
		}

	}

	return closestWaypoint;
}

int NextWaypoint(double x, double y, double theta, const vector<double> &maps_x, const vector<double> &maps_y)
{

	int closestWaypoint = ClosestWaypoint(x,y,maps_x,maps_y);

	double map_x = maps_x[closestWaypoint];
	double map_y = maps_y[closestWaypoint];

	double heading = atan2( (map_y-y),(map_x-x) );

	double angle = abs(theta-heading);

	if(angle > pi()/4)
	{
		closestWaypoint++;
	}

	return closestWaypoint;
}

// Transform from Cartesian x,y coordinates to Frenet s,d coordinates
vector<double> getFrenet(double x, double y, double theta, const vector<double> &maps_x, const vector<double> &maps_y)
{
	int next_wp = NextWaypoint(x,y, theta, maps_x,maps_y);

	int prev_wp;
	prev_wp = next_wp-1;
	if(next_wp == 0)
	{
		prev_wp  = maps_x.size()-1;
	}

	double n_x = maps_x[next_wp]-maps_x[prev_wp];
	double n_y = maps_y[next_wp]-maps_y[prev_wp];
	double x_x = x - maps_x[prev_wp];
	double x_y = y - maps_y[prev_wp];

	// find the projection of x onto n
	double proj_norm = (x_x*n_x+x_y*n_y)/(n_x*n_x+n_y*n_y);
	double proj_x = proj_norm*n_x;
	double proj_y = proj_norm*n_y;

	double frenet_d = distance(x_x,x_y,proj_x,proj_y);

	//see if d value is positive or negative by comparing it to a center point

	double center_x = 1000-maps_x[prev_wp];
	double center_y = 2000-maps_y[prev_wp];
	double centerToPos = distance(center_x,center_y,x_x,x_y);
	double centerToRef = distance(center_x,center_y,proj_x,proj_y);

	if(centerToPos <= centerToRef)
	{
		frenet_d *= -1;
	}

	// calculate s value
	double frenet_s = 0;
	for(int i = 0; i < prev_wp; i++)
	{
		frenet_s += distance(maps_x[i],maps_y[i],maps_x[i+1],maps_y[i+1]);
	}

	frenet_s += distance(0,0,proj_x,proj_y);

	return {frenet_s,frenet_d};
}

// Transform from Frenet s,d coordinates to Cartesian x,y
vector<double> getXY(double s, double d, const vector<double> &maps_s, const vector<double> &maps_x, const vector<double> &maps_y)
{
	int prev_wp = -1;

	while(s > maps_s[prev_wp+1] && (prev_wp < (int)(maps_s.size()-1) ))
	{
		prev_wp++;
	}

	int wp2 = (prev_wp+1)%maps_x.size();

	double heading = atan2((maps_y[wp2]-maps_y[prev_wp]),(maps_x[wp2]-maps_x[prev_wp]));
	// the x,y,s along the segment
	double seg_s = (s-maps_s[prev_wp]);

	double seg_x = maps_x[prev_wp]+seg_s*cos(heading);
	double seg_y = maps_y[prev_wp]+seg_s*sin(heading);

	double perp_heading = heading-pi()/2;

	double x = seg_x + d*cos(perp_heading);
	double y = seg_y + d*sin(perp_heading);

	return {x,y};
}

int main() {
  uWS::Hub h;

  // Load up map values for waypoint's x,y,s and d normalized normal vectors
  vector<double> map_waypoints_x;
  vector<double> map_waypoints_y;
  vector<double> map_waypoints_s;
  vector<double> map_waypoints_dx;
  vector<double> map_waypoints_dy;

  // Waypoint map to read from
  string map_file_ = "../data/highway_map.csv";
  // The max s value before wrapping around the track back to 0
  double max_s = 6945.554;
  
  double target_speed_mph = 0.0; //target speed at the begining [mph]
  double target_speed = target_speed_mph * 0.44704; //conversion to meters per second  
  int laneID = 1;		//include lane ID - 0, 1 or 2
  
  ifstream in_map_(map_file_.c_str(), ifstream::in);

  string line;
  while (getline(in_map_, line)) {
  	istringstream iss(line);
  	double x;
  	double y;
  	float s;
  	float d_x;
  	float d_y;
  	iss >> x;
  	iss >> y;
  	iss >> s;
  	iss >> d_x;
  	iss >> d_y;
  	map_waypoints_x.push_back(x);
  	map_waypoints_y.push_back(y);
  	map_waypoints_s.push_back(s);
  	map_waypoints_dx.push_back(d_x);
  	map_waypoints_dy.push_back(d_y);
  }
  
  h.onMessage([&laneID, &target_speed, &map_waypoints_x,&map_waypoints_y,&map_waypoints_s,&map_waypoints_dx,&map_waypoints_dy](uWS::WebSocket<uWS::SERVER> ws, char *data, size_t length, uWS::OpCode opCode) {
    // "42" at the start of the message means there's a websocket message event.
    // The 4 signifies a websocket message
    // The 2 signifies a websocket event
    //auto sdata = string(data).substr(0, length);
    //cout << sdata << endl;
    if (length && length > 2 && data[0] == '4' && data[1] == '2') {

      auto s = hasData(data);

      if (s != "") {
        auto j = json::parse(s);
        
        string event = j[0].get<string>();
        
        if (event == "telemetry") {
          // j[1] is the data JSON object
          
        	// Main car's localization Data
          	double car_x = j[1]["x"];
          	double car_y = j[1]["y"];
          	double car_s = j[1]["s"];
          	double car_d = j[1]["d"];
          	double car_yaw = j[1]["yaw"];
          	double car_speed = j[1]["speed"];

          	// Previous path data given to the Planner
          	auto previous_path_x = j[1]["previous_path_x"];
          	auto previous_path_y = j[1]["previous_path_y"];
			int previous_size = previous_path_x.size();
			
          	// Previous path's end s and d values 
          	double end_path_s = j[1]["end_path_s"];
          	double end_path_d = j[1]["end_path_d"];

          	// Sensor Fusion Data, a list of all other cars on the same side of the road.
          	auto sensor_fusion = j[1]["sensor_fusion"];
			
          	json msgJson;
			
			double number_path_points = 60;  
			double proximity_limit = 30; //meters
			double proximity_limit_behind = 10; //meters
			double speed_limit = 49;
			double lead_vehicle_speed = 49;
			int number_waypoints = 3;
			double distance_waypoints = 30;
			double lane_width = 4;
			double sim_refresh_rate = 0.02;
			int number_lanes = 3;
							
			if (previous_size > 0){
				car_s = end_path_s;
			}
				
			vector<double> next_x_vals;
          	vector<double> next_y_vals;
			
			vector<bool> freeLaneAhead;
			vector<bool> freeLaneBehind;			
						
			for (int it_lane = 0; it_lane < number_lanes; it_lane++){	
				bool it_freeLaneAhead = true;			
				bool it_freeLaneBehind = true;
				for (int it_vehicle = 0; it_vehicle < sensor_fusion.size(); it_vehicle++){
					float vehicle_d = sensor_fusion[it_vehicle][6];
					if (vehicle_d < (lane_width*(it_lane+1)) && vehicle_d > (lane_width*it_lane)){
						double vx = sensor_fusion[it_vehicle][3];
						double vy = sensor_fusion[it_vehicle][4];	
						double check_speed = sqrt(vx*vx + vy*vy);					
						double check_car_s = sensor_fusion[it_vehicle][5];
						check_car_s += (double)previous_size * sim_refresh_rate * check_speed;
						if ((check_car_s > car_s) && ((check_car_s-car_s) < proximity_limit)){
							it_freeLaneAhead = false;
							//cout << "lane " << it_lane <<" vehicle " << it_vehicle << " ahead speed: " << check_speed << endl;
							if (it_lane == laneID){
								lead_vehicle_speed = check_speed;
							}							
						}
						else if ((check_car_s < car_s) && ((car_s-check_car_s) < proximity_limit_behind)){
							it_freeLaneBehind = false;
							//cout << "lane " << it_lane <<" vehicle " << it_vehicle << " behind speed: " << check_speed << endl;
						}
					}
				}
				freeLaneAhead.push_back(it_freeLaneAhead);
				freeLaneBehind.push_back(it_freeLaneBehind);
			}
			
			bool too_close = false;
			
			//stay on the most right lane when possible & never take over on the right hand side
			//rules applied in this path planner are part of the - StVO German Autobahn driving rules
			if (laneID == 2){
				if (freeLaneAhead[2] == true){
					laneID = 2;
				}
				else if (freeLaneAhead[2] == false){
					if ((freeLaneAhead[1] == true) && (freeLaneBehind[1] == true)){
						laneID = 1;
					}
					else {
						too_close = true;
					}
				}
			}
			else if (laneID == 1){
				if (freeLaneAhead[1] == true){
					if ((freeLaneAhead[2] == true) && (freeLaneBehind[2] == true)){
						laneID = 2;
					}
					else {
						laneID = 1;
					}
				}
				else if (freeLaneAhead[1] == false){
					if ((freeLaneAhead[0] == true) && (freeLaneBehind[0] == true)) {
						laneID = 0;
					}
					else{
						too_close = true;
					}
				}
			}
			else if (laneID == 0){
				if (freeLaneAhead[0] == true) {
					if ((freeLaneAhead[1] == true) && (freeLaneBehind[1] == true)){
						laneID = 1;
					}
					else {
						laneID = 0;
					}
				}
				else if (freeLaneAhead[0] == false) {
					too_close = true;
				}
			}
			
			if (too_close){
				if (target_speed < lead_vehicle_speed){
					target_speed += 0.3;
				}
				else{
					target_speed -= 0.3;
				}
			}
			else if (target_speed < speed_limit){
				target_speed += 0.4;
			}
			
			double target_speed_met_per_sec =  target_speed*0.44704;	//conversion from mph to meters per second
			
			//equally spaced waypoints to then fit a spline to them and create a smooth path
			vector<double> points_x;
			vector<double> points_y;
			
			//reference of x, y and yaw
			double reference_x = car_x;
			double reference_y = car_y;
			double reference_yaw = deg2rad(car_yaw);
			
			if(previous_size < 2){
				double previous_car_x = car_x - cos(car_yaw);
				double previous_car_y = car_y - sin(car_yaw);
				
				points_x.push_back(previous_car_x);
				points_x.push_back(car_x);
				
				points_y.push_back(previous_car_y);
				points_y.push_back(car_y);
			}
			else {
				reference_x = previous_path_x[previous_size-1];
				reference_y = previous_path_y[previous_size-1];
				
				double reference_x_previous = previous_path_x[previous_size-2];
				double reference_y_previous = previous_path_y[previous_size-2];
				reference_yaw = atan2(reference_y-reference_y_previous, reference_x-reference_x_previous);
				
				points_x.push_back(reference_x_previous);
				points_x.push_back(reference_x);
				
				points_y.push_back(reference_y_previous);
				points_y.push_back(reference_y);					
			}
			
			
			//build few additional points to which the spline will be fit
			for(int i = 0; i < number_waypoints; i++){
				vector<double> waypoint = getXY(car_s+(i+1)*distance_waypoints, lane_width*(0.5+laneID), map_waypoints_s, map_waypoints_x, map_waypoints_y);
				points_x.push_back(waypoint[0]);
				points_y.push_back(waypoint[1]);
			}
			
			//reference frame transformation to the car coordinate system
			for(int i = 0; i < points_x.size(); i++){
				double shift_x = points_x[i] - reference_x;
				double shift_y = points_y[i] - reference_y;
				points_x[i] = (shift_x * cos(0-reference_yaw) - shift_y * sin(0-reference_yaw));
				points_y[i] = (shift_x * sin(0-reference_yaw) + shift_y * cos(0-reference_yaw));
			}
			
			spline s;
			s.set_points(points_x, points_y);		
			
			double target_x = 40.0;
			double target_y = s(target_x);
			double target_dist = sqrt(target_x*target_x+target_y*target_y);
			double x_add_on = 0;
			
			for(int i = 0; i < number_path_points; i++){
				if (i < previous_size){
					next_x_vals.push_back(previous_path_x[i]);
					next_y_vals.push_back(previous_path_y[i]);
				}
				else {
					double N = target_dist/(sim_refresh_rate*target_speed_met_per_sec);
					double x_point = x_add_on + target_x/N;
					double y_point = s(x_point);
					
					x_add_on = x_point;
					
					double x_reference = x_point;
					double y_reference = y_point;
					
					x_point = x_reference*cos(reference_yaw) - y_reference*sin(reference_yaw);
					y_point = x_reference*sin(reference_yaw) + y_reference*cos(reference_yaw);
					
					x_point += reference_x;
					y_point += reference_y;
					
					next_x_vals.push_back(x_point);
					next_y_vals.push_back(y_point);
				}
			}
			
          	msgJson["next_x"] = next_x_vals;
          	msgJson["next_y"] = next_y_vals;

          	auto msg = "42[\"control\","+ msgJson.dump()+"]";

          	//this_thread::sleep_for(chrono::milliseconds(1000));
          	ws.send(msg.data(), msg.length(), uWS::OpCode::TEXT);
          
        }
      } else {
        // Manual driving
        std::string msg = "42[\"manual\",{}]";
        ws.send(msg.data(), msg.length(), uWS::OpCode::TEXT);
      }
    }
  });

  // We don't need this since we're not using HTTP but if it's removed the
  // program doesn't compile :-(
  h.onHttpRequest([](uWS::HttpResponse *res, uWS::HttpRequest req, char *data,
                     size_t, size_t) {
    const std::string s = "<h1>Hello world!</h1>";
    if (req.getUrl().valueLength == 1) {
      res->end(s.data(), s.length());
    } else {
      // i guess this should be done more gracefully?
      res->end(nullptr, 0);
    }
  });

  h.onConnection([&h](uWS::WebSocket<uWS::SERVER> ws, uWS::HttpRequest req) {
    std::cout << "Connected!!!" << std::endl;
  });

  h.onDisconnection([&h](uWS::WebSocket<uWS::SERVER> ws, int code,
                         char *message, size_t length) {
    ws.close();
    std::cout << "Disconnected" << std::endl;
  });

  int port = 4567;
  if (h.listen(port)) {
    std::cout << "Listening to port " << port << std::endl;
  } else {
    std::cerr << "Failed to listen to port" << std::endl;
    return -1;
  }
  h.run();
}
