#pragma once

#ifndef KMEANS_CLUSTERING
#define KMEANS_CLUSTERING

#include <iostream>
#include <vector>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include <algorithm>
#include <string>

using namespace std;

template <class T>
class PointClustering
{
private:
	pair<int, int> id_point;
	int id_cluster;

	T ** trainingImage;

	int templateSizeX;
	int templateSizeY;


public:
	PointClustering(pair<int, int> id_point, T** trainingImage, int templateSizeX, int templateSizeY)
	{
		this->id_point = id_point;
		this->trainingImage = trainingImage;
		this->templateSizeX = templateSizeX;
		this->templateSizeY = templateSizeY;

		id_cluster = -1;
	}

	pair<int, int> getID()
	{
		return id_point;
	}

	int getCluster()
	{
		return id_cluster;
	}

	void setCluster(int id_cluster)
	{
		this->id_cluster = id_cluster;
	}

	T getValue(int indexX, int indexY)
	{
		return trainingImage[id_point.first + indexX][id_point.second + indexY];
	}

	void setValue(int indexX, int indexY, T value)
	{
		trainingImage[id_point.first + indexX][id_point.second + indexY] = value;
	}

	int getTotalValues()
	{
		return templateSizeX * templateSizeY;
	}

	int getSizeX()
	{
		return templateSizeX;
	}

	int getSizeY()
	{
		return templateSizeY;
	}

};

template <class T>
class ClusterClustering
{
private:
	int id_cluster;

	vector< PointClustering<T> > points;

	int sizeX;
	int sizeY;
	T** central_valuesMatrix;

public:
	ClusterClustering(int id_cluster, PointClustering<T> point)
	{
		this->id_cluster = id_cluster;
		this->sizeX = point.getSizeX();
		this->sizeY = point.getSizeY();

		central_valuesMatrix = new T*[point.getSizeX()];
		for (int i = 0; i < point.getSizeX(); i++)
		{
			central_valuesMatrix[i] = new T[point.getSizeY()];
			for (int j = 0; j < point.getSizeY(); j++)
			{
				central_valuesMatrix[i][j] = point.getValue(i, j);
			}
		}

		points.push_back(point);
	}

	void addPoint(PointClustering<T> point)
	{
		points.push_back(point);
	}

	bool removePoint(pair<int, int> id_point)
	{
		int total_points = points.size();

		for (int i = 0; i < total_points; i++)
		{
			if (points[i].getID() == id_point)				// WAIT A MINUTE!!!!!!!!! WE MUST VERIFY IT THIS LINE WORKS!!!!!!!!!!!!
			{
				points.erase(points.begin() + i);
				return true;
			}
		}
		return false;
	}

	T getCentralValue(int indexX, int indexY)
	{
		return central_valuesMatrix[indexX][indexY];
	}

	void setCentralValue(int indexX, int indexY, T value)
	{
		central_valuesMatrix[indexX][indexY] = value;
	}

	PointClustering<T> getPoint(int index)
	{
		return points[index];
	}

	int getTotalPoints()
	{
		return points.size();
	}

	int getID()
	{
		return id_cluster;
	}

	void UpdateCentroid()
	{
		for (int j = 0; j < sizeY; j++)						// for each of the Values in the Points...
		{
			for (int i = 0; i < sizeX; i++)
			{
				T sum = 0.0;									// Here we will accumulate the values of each coordinate, for all Points on the Cluster
				for (int p = 0; p < points.size(); p++)			// accumulate the values of the given 'Coordinate'
					sum += points[p].getValue(i, j);

				setCentralValue(i, j, sum / points.size());		// and put it back, averaged and updated.
			}
		}
	}

	T** getCentralValues()
	{
		return central_valuesMatrix;
	}
};

template <class T>
class KMeansClustering
{

private:
	int K; // number of clusters
	int width_TI;
	int height_TI;
	int max_iterations;		// total_points, 
	vector< ClusterClustering<T> > clusters;

	T DistanceEuclidean(ClusterClustering<T> cluster, PointClustering<T> point)
	{
		T sum = 0.0;

		for (int i = 0; i < point.getSizeX(); i++)
			for (int j = 0; j < point.getSizeY(); j++)
				sum += pow(cluster.getCentralValue(i, j) - point.getValue(i, j), 2.0);

		T leDistance = sqrt(sum);
		return leDistance;
	}

	T DistanceHamming(ClusterClustering<T> cluster, PointClustering<T> point)
	{
		T sum = 0.0;

		for (int i = 0; i < point.getSizeX(); i++)
			for (int j = 0; j < point.getSizeY(); j++)
				sum += abs(cluster.getCentralValue(i, j) - point.getValue(i, j));

		return sum;
	}

	// return ID of nearest center (uses euclidean distance)
	int getIDNearestCenter(PointClustering<T> point)
	{
		T sum = 0.0;
		T min_dist;
		int id_cluster_center = 0;
		min_dist = DistanceEuclidean(clusters[0], point);

		for (int k = 1; k < K; k++)
		{
			T dist;
			sum = 0;
			dist = DistanceEuclidean(clusters[k], point);

			if (dist < min_dist)
			{
				min_dist = dist;
				id_cluster_center = k;
			}
			if (min_dist == 0)
				break;
		}

		return id_cluster_center;
	}


public:
	KMeansClustering(int K, int max_iterations)
	{
		this->K = K;
		this->max_iterations = max_iterations;
	}

	~KMeansClustering() {
		delete this->K;
		delete this->max_iterations;
	}

	// Runs the K-means clustering algorithm over the Patterns of the given Training Image.

	vector< ClusterClustering<T> > run(T** trainingImage, int tiSizeX, int tiSizeY, int templateSizeX, int templateSizeY)
	{
		width_TI = tiSizeX - templateSizeX + 1;			// Getting the width and height to be used inside the TI.
		height_TI = tiSizeY - templateSizeY + 1;

		vector<PointClustering<T> > points = GetPointsFromTrainingImage(
			trainingImage, tiSizeX, tiSizeY, templateSizeX, templateSizeY);	// Getting all Patterns from the TI.


		vector<int> prohibited_indices;						// This structure will help us get the K initial centroids.
		for (int i = 0; i < K; i++)							// choose K distinct values for the centers of the clusters
		{
			while (true)														// to check if the coordinates weren't selected beforehand.
			{
				int id_point = rand() % points.size();							// get a Point...

				if (find(prohibited_indices.begin(), prohibited_indices.end(),	// if this point wasn't previously selected...
					id_point) == prohibited_indices.end())
				{
					{
						prohibited_indices.push_back(id_point);						// select it...
						points[id_point].setCluster(i);								// "you belong here now!"

						ClusterClustering<T> cluster(i, points[id_point]);						// and build a Cluster around it.
						clusters.push_back(cluster);								// and finally, put the newly created Cluster among his siblings.
						break; // out of the while(true) loop.

					}
				}
			}
		}


		for (int k = 0; k < K; k++)
		{
			auto leCluster = clusters[k];
			auto firstPointID = leCluster.getPoint(0).getID();
		}



		for (int iter = 1; iter <= max_iterations; iter++)							// MAIN LOOP!!! let's proceed to iterate MAX_ITERATIONS times:
		{
			cout << "--------------interation: " << iter << endl;
			bool aPointChangedCluster = false;

			for (int i = 0; i < points.size(); i++)									// FIRST PART: associates each point to the nearest CENTROID
			{
				int id_old_cluster = points[i].getCluster();							// get the previous Cluster this point was on
				int id_nearest_center = getIDNearestCenter(points[i]);					// from the updates Cluster list, get its nearest

				if (id_old_cluster != id_nearest_center)								// if they diverge...
				{
					cout << "old cluster - point[" << i << "]: " << id_old_cluster << endl;
					cout << "# new cluster: " << id_nearest_center << endl;
					if (id_old_cluster != -1)												// if the Point belonged to another Cluster...
						clusters[id_old_cluster].removePoint(points[i].getID());				// set it free!

					points[i].setCluster(id_nearest_center);								// put the Point in its nearest Cluster
					clusters[id_nearest_center].addPoint(points[i]);						// literally put it there
					aPointChangedCluster = true;											// and indicate that a Point changed its Cluster.
				}
			}




			for (int k = 0; k < K; k++)												// SECOND PART: recalculating the center of each cluster
			{
				clusters[k].UpdateCentroid();


				auto leCluster = clusters[k];
			}

			if (!aPointChangedCluster)
			{
				cout << "Break in iteration " << iter << "\n\n";
				for (int i = 0; i < points.size(); i++)									
					cout << "# Cluster from point[" << i << "]: " << points[i].getCluster() << endl;
				break;
			}
		}

		return clusters;
	}


	vector<PointClustering <T> > GetPointsFromTrainingImage(T** trainingImage, int tiSizeX, int tiSizeY, int templateSizeX, int templateSizeY)
	{
		vector<PointClustering<T> > thePoints;// (width * height);
		thePoints.reserve(width_TI * height_TI);

		for (int j = 0; j < height_TI; j++) {
			for (int i = 0; i < width_TI; i++) {
				thePoints.push_back(PointClustering<T>(make_pair(i, j), trainingImage, templateSizeX, templateSizeY)); // WHATAHELL IS HAPPENING !!!!!!!!!!!??????????????
			}
		}
		for (int i = 0; i < thePoints.size(); i++)
			cout << "Point[" << i <<"] --"<< "Pair: <" << thePoints[i].getID().first << " , " << thePoints[i].getID().second <<"> --# Values: " << thePoints[i].getValue(0, i) << " , " << thePoints[i].getValue(1, i) << " , " << thePoints[i].getValue(2, i) << " , " << thePoints[i].getValue(3, i) << endl;
		return thePoints;
	}

};

#endif // !KMEANS_CLUSTERING