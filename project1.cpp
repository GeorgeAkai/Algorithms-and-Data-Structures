//Name: George Akai
// Worked with: Riley, Elija

#include <iostream>
#include <ostream>
#include <time.h>
#include <math.h>
#include <chrono>
#include <vector>
#include <algorithm>
#include <numeric>
#include <cmath>
using namespace std;

chrono::high_resolution_clock::time_point start;
chrono::high_resolution_clock::time_point end_time;

// generate a vector of random ints in a specified range
vector<int> randomVector(int size, int low, int high)
{
  vector<int> v(size, 0);
  for (int i = 0; i < size; i++)
  {
    v[i] = rand() % (high - low + 1) + low;
    // cout << v[i] << " ";
  }
  return v;
}

// calculate sample standard deviation
double sampleSD(const vector<int> v)
{
  double mean = 0;
  for (int i = 0; i < v.size(); i++)
  {
    mean += static_cast<double>(v[i]);
  }
  mean = mean / v.size();
  double sd = 0;
  for (int i = 0; i < v.size(); i++)
  {
    sd += (static_cast<double>(v[i]) - mean) * (static_cast<double>(v[i]) - mean);
  }
  sd = sd / (v.size() - 1);
  return sqrt(sd);
}

// function to check if a vector is sorted
bool isSortedHelper(const vector<int> &v, int start, int end)
{
  if (start >= end)
  {
    // cout << "yes" << endl;
    return true;
  }
  if (v[start] > v[end])
  {
    // cout << "nope" << endl;
    return false;
  }
  return isSortedHelper(v, start + 1, end);
}

bool isSorted(const vector<int> &v)
{
  return isSortedHelper(v, 0, v.size() - 1);
}

// Insertion Sort Implementation
void InsertionSort(vector<int> &v)
{
  int i = 1;
  while (i < v.size())
  {
    int j = i;
    while (j > 0 && v[j] < v[j - 1])
    {
      int temp = v[j - 1];
      v[j - 1] = v[j];
      v[j] = temp;
      j--;
    }
    i++;
  }
}

// implementing selection sort
void SelectionSort(vector<int> &v)
{
  unsigned size = v.size();
  for (unsigned int i = 0; i <= size; i++)
  {
    unsigned int small = i;
    for (unsigned int j = i + 1; j <= size - 1; j++)
    {
      if (v[j] < v[small])
      {
        small = j;
      }
    }
    if (small != i)
    {
      int temp = v[i];
      v[i] = v[small];
      v[small] = temp;
    }
  }
}
// implementing bubble sort
void BubbleSort(vector<int> &v)
{
  bool sorted = false;
  int size = v.size();
  while (!sorted)
  {
    sorted = true;
    for (int i = 1; i < size; i++)
    {
      if (v[i - 1] > v[i])
      {
        int temp = v[i - 1];
        v[i - 1] = v[i];
        v[i] = temp;
        sorted = false;
      }
    }
    size--;
  }
}

// implementing quick sort
void QuickSort(vector<int> &v)
{
  if (v.size() <= 1)
  {
    return;
  }
  int pivot = v[0];
  vector<int> low;
  vector<int> high;
  for (int i = 1; i < v.size(); i++)
  {
    if (v[i] <= pivot)
    {
      low.push_back(v[i]);
    }
    else
    {
      high.push_back(v[i]);
    }
  }
  QuickSort(low);
  QuickSort(high);

  int index = 0;

  for (int i = 0; i < low.size(); i++)
  {
    v[index++] = low[i];
  }

  v[index++] = pivot;

  for (int i = 0; i < high.size(); i++)
  {
    v[index++] = high[i];
  }
}

// Reversely sort through a vector using bubblesort
void ReverseBubbleSort(vector<int> &v)
{
  bool sorted = false;
  int size = v.size();
  while (!sorted)
  {
    sorted = true;
    for (int i = 1; i < size; i++)
    {
      if (v[i - 1] < v[i])
      {
        int temp = v[i - 1];
        v[i - 1] = v[i];
        v[i] = temp;
        sorted = false;
      }
    }
    size--;
  }
}

// Generating  a function for the Best Case for Quick Sort
std::vector<int> BalanceOrder(std::vector<int> &v)
{
  if (v.size() <= 1)
  {
    return v;
  }
  int mid = v.size() / 2; // choose the median as the pivot

  std::vector<int> output = {v[mid]};
  std::vector<int> left(v.begin(), v.begin() + mid);    // rec call for left subarr
  std::vector<int> right(v.begin() + mid + 1, v.end()); // recursive call for right subarray

  // add them together
  std::vector<int> balanceleft = BalanceOrder(left);
  std::vector<int> balanceright = BalanceOrder(right);

  // append to the output
  output.insert(output.end(), balanceleft.begin(), balanceleft.end());
  output.insert(output.end(), balanceright.begin(), balanceright.end());

  // return the final output
  return output;
}

int main()
{
  // Testing and Experiments: Question 1
  cout << "QUESTION 1" << endl;
  srand(time(NULL));
  // timing functions in c++
  // store the runtimes
  vector<double> runtime;
  // Run Sorting Algorithms on 10 random vectors of size 100, collects the runtimes, successful sorting, print min, mean, sd, max
  for (int i = 0; i < 10; i++)
  {
    vector<int> v = randomVector(100, 0, 99);
    start = chrono::high_resolution_clock::now();
    BubbleSort(v);
    end_time = chrono::high_resolution_clock::now();
    double time = chrono::duration_cast<chrono::duration<double>>(end_time - start).count();
    runtime.push_back(time);

    // calc min runtime
    double min = *min_element(runtime.begin(), runtime.end());

    // calc the mean runtime
    double mean = accumulate(runtime.begin(), runtime.end(), 0.0);

    // calc the maximum runtime
    double max = *max_element(runtime.begin(), runtime.end());

    // output the results of Bubble Sort
    cout << "*******************" << endl;
    cout << "Bubble Sort on 10 Vectors of length 100" << endl;
    if (isSorted(v))
    {
      cout << "Sorting successful!" << endl;
    }
    cout << "Minimum: " << min << " sec; " << "Mean: " << mean << " sec; " << "Standard Deviation: " << sampleSD(v) << " sec; " << "Maximum: " << max << " sec" << endl;
    cout << "*******************" << endl;
  }

  // Selection Sort
  for (int i = 0; i < 10; i++)
  {
    vector<int> v = randomVector(100, 0, 99);
    start = chrono::high_resolution_clock::now();
    SelectionSort(v);
    end_time = chrono::high_resolution_clock::now();
    double time = chrono::duration_cast<chrono::duration<double>>(end_time - start).count();
    runtime.push_back(time);

    // calc min runtime
    double min = *min_element(runtime.begin(), runtime.end());

    // calc the mean runtime
    double mean = accumulate(runtime.begin(), runtime.end(), 0.0);

    // calc the maximum runtime
    double max = *max_element(runtime.begin(), runtime.end());

    // output the results of Selection Sort
    cout << "*******************" << endl;
    cout << "Selection Sort on 10 Vectors of length 100" << endl;
    if (isSorted(v))
    {
      cout << "Sorting successful!" << endl;
    }
    cout << "Minimum: " << min << " sec; " << "Mean: " << mean << " sec; " << "Standard Deviation: " << sampleSD(v) << " sec; " << "Maximum: " << max << " sec" << endl;
    cout << "*******************" << endl;
  }

  // Insertion sort runtimes
  for (int i = 0; i < 10; i++)
  {
    vector<int> v = randomVector(100, 0, 99);
    start = chrono::high_resolution_clock::now();
    InsertionSort(v);
    end_time = chrono::high_resolution_clock::now();
    double time = chrono::duration_cast<chrono::duration<double>>(end_time - start).count();
    runtime.push_back(time);

    // calc min runtime
    double min = *min_element(runtime.begin(), runtime.end());

    // calc the mean runtime
    double mean = accumulate(runtime.begin(), runtime.end(), 0.0);

    // calc the maximum runtime
    double max = *max_element(runtime.begin(), runtime.end());

    // output the results of Insertion Sort
    cout << "*******************" << endl;
    cout << "Insertion Sort on 10 Vectors of length 100" << endl;
    if (isSorted(v))
    {
      cout << "Sorting successful!" << endl;
    }
    cout << "Minimum: " << min << " sec; " << "Mean: " << mean << " sec; " << "Standard Deviation: " << sampleSD(v) << " sec; " << "Maximum: " << max << " sec" << endl;
    cout << "*******************" << endl;
  }

  // Quick Sort runtimes
  for (int i = 0; i < 10; i++)
  {
    vector<int> v = randomVector(100, 0, 99);
    start = chrono::high_resolution_clock::now();
    QuickSort(v);
    end_time = chrono::high_resolution_clock::now();
    double time = chrono::duration_cast<chrono::duration<double>>(end_time - start).count();
    runtime.push_back(time);

    // calc min runtime
    double min = *min_element(runtime.begin(), runtime.end());

    // calc the mean runtime
    double mean = accumulate(runtime.begin(), runtime.end(), 0.0);

    // calc the maximum runtime
    double max = *max_element(runtime.begin(), runtime.end());

    // output the results of Quick Sort
    cout << "*******************" << endl;
    cout << "Quick Sort on 10 Vectors of length 100" << endl;
    if (isSorted(v))
    {
      cout << "Sorting successful!" << endl;
    }
    cout << "Minimum: " << min << " sec; " << "Mean: " << mean << " sec; " << "Standard Deviation: " << sampleSD(v) << " sec; " << "Maximum: " << max << " sec" << endl;
    cout << "*******************" << endl;
  }

  // Testing and Experiments: Question 2
  cout << "QUESTION 2" << endl;
  // Generating average case vectors
  // generating 50 random average case vectors of length 10
  cout << "Avg case for size 10 vector" << endl;
  for (int i = 0; i < 50; i++)
  {
    vector<int> v = randomVector(10, 0, 99);
    start = chrono::high_resolution_clock::now();
    BubbleSort(v);
    end_time = chrono::high_resolution_clock::now();
    double time = chrono::duration_cast<chrono::duration<double>>(end_time - start).count();
    cout << "Bubble Sort average runtime (size 10): " << time << endl;
  }
  cout << endl;

  for (int i = 0; i < 50; i++)
  {
    vector<int> v = randomVector(10, 0, 99);
    start = chrono::high_resolution_clock::now();
    SelectionSort(v);
    end_time = chrono::high_resolution_clock::now();
    double time = chrono::duration_cast<chrono::duration<double>>(end_time - start).count();
    cout << "Selection Sort average runtime (size 10): " << time << endl;
  }
  cout << endl;

  for (int i = 0; i < 50; i++)
  {
    vector<int> v = randomVector(10, 0, 99);
    start = chrono::high_resolution_clock::now();
    InsertionSort(v);
    end_time = chrono::high_resolution_clock::now();
    double time = chrono::duration_cast<chrono::duration<double>>(end_time - start).count();
    cout << "Insertion Sort average runtime (size 10): " << time << endl;
  }
  cout << endl;

  for (int i = 0; i < 50; i++)
  {
    vector<int> v = randomVector(10, 0, 99);
    start = chrono::high_resolution_clock::now();
    QuickSort(v);
    end_time = chrono::high_resolution_clock::now();
    double time = chrono::duration_cast<chrono::duration<double>>(end_time - start).count();
    cout << "Quick Sort average runtime (size 10): " << time << endl;
  }
  cout << endl;

  // Generating best-case runtimes for each algorithm
  cout << "Best case for length 10" << endl;
  for (int i = 0; i < 50; i++)
  {
    vector<int> v = randomVector(10, 0, 99);
    // Sort the vector for best case runtime.
    SelectionSort(v);

    // bubble sort best case runtime
    start = chrono::high_resolution_clock::now();
    BubbleSort(v);
    end_time = chrono::high_resolution_clock::now();
    double time = chrono::duration_cast<chrono::duration<double>>(end_time - start).count();
    cout << "Bubble sort best-case runtime (size 10): " << time << endl;
    // selection sort
    start = chrono::high_resolution_clock::now();
    SelectionSort(v);
    end_time = chrono::high_resolution_clock::now();
    double ssTime = chrono::duration_cast<chrono::duration<double>>(end_time - start).count();
    cout << "Selection sort best-case runtime (size 10): " << ssTime << endl;
    // insertion sort
    start = chrono::high_resolution_clock::now();
    InsertionSort(v);
    end_time = chrono::high_resolution_clock::now();
    double isTime = chrono::duration_cast<chrono::duration<double>>(end_time - start).count();
    cout << "Insertion sort best-case runtime (size 10): " << isTime << endl;
    // Quick sort
    start = chrono::high_resolution_clock::now();
    QuickSort(v);
    end_time = chrono::high_resolution_clock::now();
    double qsTime = chrono::duration_cast<chrono::duration<double>>(end_time - start).count();
    cout << "Quick sort worst-case runtime (size 10): " << qsTime << endl;
  }
  cout << endl;

  // Reversely sorted Vectors for worst runtimes for bb, ins, sel
  cout << "Worst case for 10" << endl;
  for (int i = 0; i < 50; i++)
  {
    vector<int> v = randomVector(10, 0, 99);
    // reversely sort the vector for worst case runtime.
    ReverseBubbleSort(v);

    // bubble sort worst-case runtime
    start = chrono::high_resolution_clock::now();
    BubbleSort(v);
    end_time = chrono::high_resolution_clock::now();
    double time = chrono::duration_cast<chrono::duration<double>>(end_time - start).count();
    cout << "Bubble sort worst-case runtime (size 10): " << time << endl;
    // selection sort
    start = chrono::high_resolution_clock::now();
    SelectionSort(v);
    end_time = chrono::high_resolution_clock::now();
    double ssTime = chrono::duration_cast<chrono::duration<double>>(end_time - start).count();
    cout << "Selection sort worst-case runtime (size 10): " << ssTime << endl;
    // insertion sort
    start = chrono::high_resolution_clock::now();
    InsertionSort(v);
    end_time = chrono::high_resolution_clock::now();
    double isTime = chrono::duration_cast<chrono::duration<double>>(end_time - start).count();
    cout << "Insertion sort worst-case runtime (size 10): " << isTime << endl;
    // Quick sort
    start = chrono::high_resolution_clock::now();
    QuickSort(v);
    end_time = chrono::high_resolution_clock::now();
    double qsTime = chrono::duration_cast<chrono::duration<double>>(end_time - start).count();
    cout << "Quick sort worst-case runtime (size 10): " << qsTime << endl;
  }
  cout << endl;

  // Repeat for size 100
  // generating 50 random average case vectors of length 10
  cout << "Avg case for 100" << endl;
  for (int i = 0; i < 50; i++)
  {
    vector<int> v = randomVector(100, 0, 99);
    start = chrono::high_resolution_clock::now();
    BubbleSort(v);
    end_time = chrono::high_resolution_clock::now();
    double time = chrono::duration_cast<chrono::duration<double>>(end_time - start).count();
    cout << "Bubble Sort average runtime (size 100): " << time << endl;
  }
  cout << endl;

  for (int i = 0; i < 50; i++)
  {
    vector<int> v = randomVector(100, 0, 99);
    start = chrono::high_resolution_clock::now();
    SelectionSort(v);
    end_time = chrono::high_resolution_clock::now();
    double time = chrono::duration_cast<chrono::duration<double>>(end_time - start).count();
    cout << "Selection Sort average runtime (size 100): " << time << endl;
  }
  cout << endl;

  for (int i = 0; i < 50; i++)
  {
    vector<int> v = randomVector(100, 0, 99);
    start = chrono::high_resolution_clock::now();
    InsertionSort(v);
    end_time = chrono::high_resolution_clock::now();
    double time = chrono::duration_cast<chrono::duration<double>>(end_time - start).count();
    cout << "Insertion Sort average runtime (size 100): " << time << endl;
  }
  cout << endl;

  for (int i = 0; i < 50; i++)
  {
    vector<int> v = randomVector(100, 0, 99);
    start = chrono::high_resolution_clock::now();
    QuickSort(v);
    end_time = chrono::high_resolution_clock::now();
    double time = chrono::duration_cast<chrono::duration<double>>(end_time - start).count();
    cout << "Quick Sort average runtime (size 100): " << time << endl;
  }
  cout << endl;

  // Generating best-case runtimes for each algorithm
  cout << "Best case for 100" << endl;
  for (int i = 0; i < 50; i++)
  {
    vector<int> v = randomVector(100, 0, 99);
    // Sort the vector for best case runtime.
    SelectionSort(v);

    // bubble sort best case runtime
    start = chrono::high_resolution_clock::now();
    BubbleSort(v);
    end_time = chrono::high_resolution_clock::now();
    double time = chrono::duration_cast<chrono::duration<double>>(end_time - start).count();
    cout << "Bubble sort best-case runtime (size 100): " << time << endl;
    // selection sort
    start = chrono::high_resolution_clock::now();
    SelectionSort(v);
    end_time = chrono::high_resolution_clock::now();
    double ssTime = chrono::duration_cast<chrono::duration<double>>(end_time - start).count();
    cout << "Selection sort best-case runtime (size 100): " << ssTime << endl;
    // insertion sort
    start = chrono::high_resolution_clock::now();
    InsertionSort(v);
    end_time = chrono::high_resolution_clock::now();
    double isTime = chrono::duration_cast<chrono::duration<double>>(end_time - start).count();
    cout << "Insertion sort best-case runtime (size 100): " << isTime << endl;
    // Quick sort
    start = chrono::high_resolution_clock::now();
    QuickSort(v);
    end_time = chrono::high_resolution_clock::now();
    double qsTime = chrono::duration_cast<chrono::duration<double>>(end_time - start).count();
    cout << "Quick sort worst-case runtime (size 100): " << qsTime << endl;
  }
  cout << endl;

  // Reversely sorted Vectors for worst runtimes for bb, ins, sel
  cout << "Worst case for 100" << endl;
  for (int i = 0; i < 50; i++)
  {
    vector<int> v = randomVector(100, 0, 99);
    // reversely sort the vector for worst case runtime.
    ReverseBubbleSort(v);

    // bubble sort worst-case runtime
    start = chrono::high_resolution_clock::now();
    BubbleSort(v);
    end_time = chrono::high_resolution_clock::now();
    double time = chrono::duration_cast<chrono::duration<double>>(end_time - start).count();
    cout << "Bubble sort worst-case runtime (size 100): " << time << endl;
    // selection sort
    start = chrono::high_resolution_clock::now();
    SelectionSort(v);
    end_time = chrono::high_resolution_clock::now();
    double ssTime = chrono::duration_cast<chrono::duration<double>>(end_time - start).count();
    cout << "Selection sort worst-case runtime (size 100): " << ssTime << endl;
    // insertion sort
    start = chrono::high_resolution_clock::now();
    InsertionSort(v);
    end_time = chrono::high_resolution_clock::now();
    double isTime = chrono::duration_cast<chrono::duration<double>>(end_time - start).count();
    cout << "Insertion sort worst-case runtime (size 100): " << isTime << endl;
    // Quick sort
    start = chrono::high_resolution_clock::now();
    QuickSort(v);
    end_time = chrono::high_resolution_clock::now();
    double qsTime = chrono::duration_cast<chrono::duration<double>>(end_time - start).count();
    cout << "Quick sort worst-case runtime (size 100): " << qsTime << endl;
  }
  cout << endl;

  // size = 1000
  // generating 50 random average case vectors of length 10
  cout << "Avg case for 1000" << endl;
  for (int i = 0; i < 50; i++)
  {
    vector<int> v = randomVector(1000, 0, 999);
    start = chrono::high_resolution_clock::now();
    BubbleSort(v);
    end_time = chrono::high_resolution_clock::now();
    double time = chrono::duration_cast<chrono::duration<double>>(end_time - start).count();
    cout << "Bubble Sort average runtime (size 1000): " << time << endl;
  }
  cout << endl;

  for (int i = 0; i < 50; i++)
  {
    vector<int> v = randomVector(1000, 0, 999);
    start = chrono::high_resolution_clock::now();
    SelectionSort(v);
    end_time = chrono::high_resolution_clock::now();
    double time = chrono::duration_cast<chrono::duration<double>>(end_time - start).count();
    cout << "Selection Sort average runtime (size 1000): " << time << endl;
  }
  cout << endl;

  for (int i = 0; i < 50; i++)
  {
    vector<int> v = randomVector(1000, 0, 999);
    start = chrono::high_resolution_clock::now();
    InsertionSort(v);
    end_time = chrono::high_resolution_clock::now();
    double time = chrono::duration_cast<chrono::duration<double>>(end_time - start).count();
    cout << "Insertion Sort average runtime (size 1000): " << time << endl;
  }
  cout << endl;

  for (int i = 0; i < 50; i++)
  {
    vector<int> v = randomVector(1000, 0, 999);
    start = chrono::high_resolution_clock::now();
    QuickSort(v);
    end_time = chrono::high_resolution_clock::now();
    double time = chrono::duration_cast<chrono::duration<double>>(end_time - start).count();
    cout << "Quick Sort average runtime (size 1000): " << time << endl;
  }
  cout << endl;

  // Generating best-case runtimes for each algorithm
  cout << "Best case for 1000" << endl;
  for (int i = 0; i < 50; i++)
  {
    vector<int> v = randomVector(1000, 0, 999);
    // Sort the vector for best case runtime.
    SelectionSort(v);

    // bubble sort best case runtime
    start = chrono::high_resolution_clock::now();
    BubbleSort(v);
    end_time = chrono::high_resolution_clock::now();
    double time = chrono::duration_cast<chrono::duration<double>>(end_time - start).count();
    cout << "Bubble sort best-case runtime (size 1000): " << time << endl;
    // selection sort
    start = chrono::high_resolution_clock::now();
    SelectionSort(v);
    end_time = chrono::high_resolution_clock::now();
    double ssTime = chrono::duration_cast<chrono::duration<double>>(end_time - start).count();
    cout << "Selection sort best-case runtime (size 1000): " << ssTime << endl;
    // insertion sort
    start = chrono::high_resolution_clock::now();
    InsertionSort(v);
    end_time = chrono::high_resolution_clock::now();
    double isTime = chrono::duration_cast<chrono::duration<double>>(end_time - start).count();
    cout << "Insertion sort best-case runtime (size 1000): " << isTime << endl;
    // Quick sort
    start = chrono::high_resolution_clock::now();
    QuickSort(v);
    end_time = chrono::high_resolution_clock::now();
    double qsTime = chrono::duration_cast<chrono::duration<double>>(end_time - start).count();
    cout << "Quick sort worst-case runtime (size 1000): " << qsTime << endl;
  }
  cout << endl;

  // Reversely sorted Vectors for worst runtimes for bb, ins, sel
  cout << "Worst case for 1000" << endl;
  for (int i = 0; i < 50; i++)
  {
    vector<int> v = randomVector(1000, 0, 999);
    // reversely sort the vector for worst case runtime.
    ReverseBubbleSort(v);

    // bubble sort worst-case runtime
    start = chrono::high_resolution_clock::now();
    BubbleSort(v);
    end_time = chrono::high_resolution_clock::now();
    double time = chrono::duration_cast<chrono::duration<double>>(end_time - start).count();
    cout << "Bubble sort worst-case runtime (size 1000): " << time << endl;
    // selection sort
    start = chrono::high_resolution_clock::now();
    SelectionSort(v);
    end_time = chrono::high_resolution_clock::now();
    double ssTime = chrono::duration_cast<chrono::duration<double>>(end_time - start).count();
    cout << "Selection sort worst-case runtime (size 1000): " << ssTime << endl;
    // insertion sort
    start = chrono::high_resolution_clock::now();
    InsertionSort(v);
    end_time = chrono::high_resolution_clock::now();
    double isTime = chrono::duration_cast<chrono::duration<double>>(end_time - start).count();
    cout << "Insertion sort worst-case runtime (size 1000): " << isTime << endl;
    // Quick sort
    start = chrono::high_resolution_clock::now();
    QuickSort(v);
    end_time = chrono::high_resolution_clock::now();
    double qsTime = chrono::duration_cast<chrono::duration<double>>(end_time - start).count();
    cout << "Quick sort worst-case runtime (size 1000): " << qsTime << endl;
  }
  cout << endl;

  // size of Vector  = 5000
  cout << "Average case for 5000" << endl;
  for (int i = 0; i < 50; i++)
  {
    vector<int> v = randomVector(5000, 0, 4999);
    start = chrono::high_resolution_clock::now();
    BubbleSort(v);
    end_time = chrono::high_resolution_clock::now();
    double time = chrono::duration_cast<chrono::duration<double>>(end_time - start).count();
    cout << "Bubble Sort average runtime (size 5000): " << time << endl;
  }
  cout << endl;

  for (int i = 0; i < 50; i++)
  {
    vector<int> v = randomVector(5000, 0, 4999);
    start = chrono::high_resolution_clock::now();
    SelectionSort(v);
    end_time = chrono::high_resolution_clock::now();
    double time = chrono::duration_cast<chrono::duration<double>>(end_time - start).count();
    cout << "Selection Sort average runtime (size 5000): " << time << endl;
  }
  cout << endl;

  for (int i = 0; i < 50; i++)
  {
    vector<int> v = randomVector(5000, 0, 4999);
    start = chrono::high_resolution_clock::now();
    InsertionSort(v);
    end_time = chrono::high_resolution_clock::now();
    double time = chrono::duration_cast<chrono::duration<double>>(end_time - start).count();
    cout << "Insertion Sort average runtime (size 1000): " << time << endl;
  }
  cout << endl;

  for (int i = 0; i < 50; i++)
  {
    vector<int> v = randomVector(5000, 0, 4999);
    start = chrono::high_resolution_clock::now();
    QuickSort(v);
    end_time = chrono::high_resolution_clock::now();
    double time = chrono::duration_cast<chrono::duration<double>>(end_time - start).count();
    cout << "Quick Sort average runtime (size 5000): " << time << endl;
  }
  cout << endl;

  // Generating best-case runtimes for each algorithm
  cout << "Best Case for 5000 elements" << endl;
  for (int i = 0; i < 50; i++)
  {
    vector<int> v = randomVector(5000, 0, 4999);
    // Sort the vector for best case runtime.
    SelectionSort(v);

    // bubble sort best case runtime
    start = chrono::high_resolution_clock::now();
    BubbleSort(v);
    end_time = chrono::high_resolution_clock::now();
    double time = chrono::duration_cast<chrono::duration<double>>(end_time - start).count();
    cout << "Bubble sort best-case runtime (size 5000): " << time << endl;
    // selection sort
    start = chrono::high_resolution_clock::now();
    SelectionSort(v);
    end_time = chrono::high_resolution_clock::now();
    double ssTime = chrono::duration_cast<chrono::duration<double>>(end_time - start).count();
    cout << "Selection sort best-case runtime (size 5000): " << ssTime << endl;
    // insertion sort
    start = chrono::high_resolution_clock::now();
    InsertionSort(v);
    end_time = chrono::high_resolution_clock::now();
    double isTime = chrono::duration_cast<chrono::duration<double>>(end_time - start).count();
    cout << "Insertion sort best-case runtime (size 5000): " << isTime << endl;
    // Quick sort
    start = chrono::high_resolution_clock::now();
    QuickSort(v);
    end_time = chrono::high_resolution_clock::now();
    double qsTime = chrono::duration_cast<chrono::duration<double>>(end_time - start).count();
    cout << "Quick sort worst-case runtime (size 5000): " << qsTime << endl;
  }
  cout << "Worst Cases" << endl;
  // Reversely sorted Vectors for worst runtimes for bb, ins, sel
  for (int i = 0; i < 50; i++)
  {
    vector<int> v = randomVector(5000, 0, 4999);
    // reversely sort the vector for worst case runtime.
    ReverseBubbleSort(v);

    // bubble sort worst-case runtime
    start = chrono::high_resolution_clock::now();
    BubbleSort(v);
    end_time = chrono::high_resolution_clock::now();
    double time = chrono::duration_cast<chrono::duration<double>>(end_time - start).count();
    cout << "Bubble sort worst-case runtime (size 5000): " << time << endl;
    // selection sort
    start = chrono::high_resolution_clock::now();
    SelectionSort(v);
    end_time = chrono::high_resolution_clock::now();
    double ssTime = chrono::duration_cast<chrono::duration<double>>(end_time - start).count();
    cout << "Selection sort worst-case runtime (size 5000): " << ssTime << endl;
    // insertion sort
    start = chrono::high_resolution_clock::now();
    InsertionSort(v);
    end_time = chrono::high_resolution_clock::now();
    double isTime = chrono::duration_cast<chrono::duration<double>>(end_time - start).count();
    cout << "Insertion sort worst-case runtime (size 5000): " << isTime << endl;
    // Quick sort
    start = chrono::high_resolution_clock::now();
    QuickSort(v);
    end_time = chrono::high_resolution_clock::now();
    double qsTime = chrono::duration_cast<chrono::duration<double>>(end_time - start).count();
    cout << "Quick sort worst-case runtime (size 5000): " << qsTime << endl;
  }
  cout << endl;

  // size of 10000 elements in a vector
  cout << "Average Case Runtime for 10000 elements" << endl;
  for (int i = 0; i < 50; i++)
  {
    vector<int> v = randomVector(10000, 0, 9999);
    start = chrono::high_resolution_clock::now();
    BubbleSort(v);
    end_time = chrono::high_resolution_clock::now();
    double time = chrono::duration_cast<chrono::duration<double>>(end_time - start).count();
    cout << "Bubble Sort average runtime (size 10000): " << time << endl;
  }
  cout << endl;

  for (int i = 0; i < 50; i++)
  {
    vector<int> v = randomVector(10000, 0, 9999);
    start = chrono::high_resolution_clock::now();
    SelectionSort(v);
    end_time = chrono::high_resolution_clock::now();
    double time = chrono::duration_cast<chrono::duration<double>>(end_time - start).count();
    cout << "Selection Sort average runtime (size 10000): " << time << endl;
  }
  cout << endl;

  for (int i = 0; i < 50; i++)
  {
    vector<int> v = randomVector(10000, 0, 9999);
    start = chrono::high_resolution_clock::now();
    InsertionSort(v);
    end_time = chrono::high_resolution_clock::now();
    double time = chrono::duration_cast<chrono::duration<double>>(end_time - start).count();
    cout << "Insertion Sort average runtime (size 10000): " << time << endl;
  }
  cout << endl;

  for (int i = 0; i < 50; i++)
  {
    vector<int> v = randomVector(10000, 0, 9999);
    start = chrono::high_resolution_clock::now();
    QuickSort(v);
    end_time = chrono::high_resolution_clock::now();
    double time = chrono::duration_cast<chrono::duration<double>>(end_time - start).count();
    cout << "Quick Sort average runtime (size 10000): " << time << endl;
  }
  cout << endl;

  // Generating best-case runtimes for each algorithm
  cout << "Best Cases" << endl;
  for (int i = 0; i < 50; i++)
  {
    vector<int> v = randomVector(10000, 0, 9999);
    // Sort the vector for best case runtime.
    SelectionSort(v);

    // bubble sort best case runtime
    start = chrono::high_resolution_clock::now();
    BubbleSort(v);
    end_time = chrono::high_resolution_clock::now();
    double time = chrono::duration_cast<chrono::duration<double>>(end_time - start).count();
    cout << "Bubble sort best-case runtime (size 10000): " << time << endl;
    // selection sort
    start = chrono::high_resolution_clock::now();
    SelectionSort(v);
    end_time = chrono::high_resolution_clock::now();
    double ssTime = chrono::duration_cast<chrono::duration<double>>(end_time - start).count();
    cout << "Selection sort best-case runtime (size 10000): " << ssTime << endl;
    // insertion sort
    start = chrono::high_resolution_clock::now();
    InsertionSort(v);
    end_time = chrono::high_resolution_clock::now();
    double isTime = chrono::duration_cast<chrono::duration<double>>(end_time - start).count();
    cout << "Insertion sort best-case runtime (size 10000): " << isTime << endl;
    // Quick sort
    start = chrono::high_resolution_clock::now();
    QuickSort(v);
    end_time = chrono::high_resolution_clock::now();
    double qsTime = chrono::duration_cast<chrono::duration<double>>(end_time - start).count();
    cout << "Quick sort worst-case runtime (size 10000): " << qsTime << endl;
  }
  cout << "Worst Cases" << endl;
  // Reversely sorted Vectors for worst runtimes for bb, ins, sel
  for (int i = 0; i < 50; i++)
  {
    vector<int> v = randomVector(10000, 0, 9999);
    // reversely sort the vector for worst case runtime.
    ReverseBubbleSort(v);

    // bubble sort worst-case runtime
    start = chrono::high_resolution_clock::now();
    BubbleSort(v);
    end_time = chrono::high_resolution_clock::now();
    double time = chrono::duration_cast<chrono::duration<double>>(end_time - start).count();
    cout << "Bubble sort worst-case runtime (size 10000): " << time << endl;
    // selection sort
    start = chrono::high_resolution_clock::now();
    SelectionSort(v);
    end_time = chrono::high_resolution_clock::now();
    double ssTime = chrono::duration_cast<chrono::duration<double>>(end_time - start).count();
    cout << "Selection sort worst-case runtime (size 10000): " << ssTime << endl;
    // insertion sort
    start = chrono::high_resolution_clock::now();
    InsertionSort(v);
    end_time = chrono::high_resolution_clock::now();
    double isTime = chrono::duration_cast<chrono::duration<double>>(end_time - start).count();
    cout << "Insertion sort worst-case runtime (size 10000): " << isTime << endl;
    // Quick sort
    start = chrono::high_resolution_clock::now();
    QuickSort(v);
    end_time = chrono::high_resolution_clock::now();
    double qsTime = chrono::duration_cast<chrono::duration<double>>(end_time - start).count();
    cout << "Quick sort worst-case runtime (size 10000): " << qsTime << endl;
  }
  cout << endl;

  // Best case for QuickSort
  cout << "Best Case for QuickSort" << endl;
  for (int i = 0; i < 50; i++)
  {
    vector<int> v = randomVector(10, 0, 9);
    BubbleSort(v);       // I am just sorting the vector here using bubblesort
    v = BalanceOrder(v); // balance the vector for quicksort best case
    start = chrono::high_resolution_clock::now();
    QuickSort(v); // do quicksort in the balanced vector
    end_time = chrono::high_resolution_clock::now();
    double bestTime = chrono::duration_cast<chrono::duration<double>>(end_time - start).count();
    cout << "Quick sort BEST-case runtime (size 10): " << bestTime << endl;
  }
  cout << endl;

  for (int i = 0; i < 50; i++)
  {
    vector<int> v = randomVector(100, 0, 99);
    BubbleSort(v);       // I am just sorting the vector here using bubblesort
    v = BalanceOrder(v); // balance the vector for quicksort best case
    start = chrono::high_resolution_clock::now();
    QuickSort(v); // do quicksort in the balanced vector
    end_time = chrono::high_resolution_clock::now();
    double bestTime = chrono::duration_cast<chrono::duration<double>>(end_time - start).count();
    cout << "Quick sort BEST-case runtime (size 100): " << bestTime << endl;
  }
  cout << endl;
  for (int i = 0; i < 50; i++)
  {
    vector<int> v = randomVector(1000, 0, 999);
    BubbleSort(v);       // I am just sorting the vector here using bubblesort
    v = BalanceOrder(v); // balance the vector for quicksort best case
    start = chrono::high_resolution_clock::now();
    QuickSort(v); // do quicksort in the balanced vector
    end_time = chrono::high_resolution_clock::now();
    double bestTime = chrono::duration_cast<chrono::duration<double>>(end_time - start).count();
    cout << "Quick sort BEST-case runtime (size 1000): " << bestTime << endl;
  }
  cout << endl;
  for (int i = 0; i < 50; i++)
  {
    vector<int> v = randomVector(5000, 0, 4999);
    BubbleSort(v);       // I am just sorting the vector here using bubblesort
    v = BalanceOrder(v); // balance the vector for quicksort best case
    start = chrono::high_resolution_clock::now();
    QuickSort(v); // do quicksort in the balanced vector
    end_time = chrono::high_resolution_clock::now();
    double bestTime = chrono::duration_cast<chrono::duration<double>>(end_time - start).count();
    cout << "Quick sort BEST-case runtime (size 5000): " << bestTime << endl;
  }
  cout << endl;
  for (int i = 0; i < 50; i++)
  {
    vector<int> v = randomVector(10000, 0, 9999);
    BubbleSort(v);       // I am just sorting the vector here using bubblesort
    v = BalanceOrder(v); // balance the vector for quicksort best case
    start = chrono::high_resolution_clock::now();
    QuickSort(v); // do quicksort in the balanced vector
    end_time = chrono::high_resolution_clock::now();
    double bestTime = chrono::duration_cast<chrono::duration<double>>(end_time - start).count();
    cout << "Quick sort BEST-case runtime (size 10000): " << bestTime << endl;
  }
  cout << endl;

  // END OF THE PROGRAM
  return 0;
}