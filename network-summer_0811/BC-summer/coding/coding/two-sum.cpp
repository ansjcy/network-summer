////
////  main.cpp
////  coding
////
////  Created by Anakin on 16/7/20.
////  Copyright © 2016年 Anakin. All rights reserved.
////
//
//#include <iostream>
//#include <math.h>
//#include <vector>
//#include <map>
//using namespace std;
//
//class Solution {
//public:
////    vector<int> twoSum(vector<int>& nums, int target) {
////        vector<int> result;
////        
////        for(int i = 0; i < nums.size(); i++)
////        {
////            for(int j = i+1; j < nums.size(); j++)
////            {
////                if(nums[i] + nums[j] == target)
////                {
////                    result.push_back(i);
////                    result.push_back(j);
////                    return result;
////                }
////            }
////        }
////        return result;
////    }
//    vector<int> twoSum(vector<int>& nums, int target) {
//        vector<int> result;
//        map<int, int> saved;
//        for(int i = 0; i < nums.size(); i++)
//        {
//            if(saved.find(target-nums[i]) != saved.end())
//            {
//                result.push_back(saved[target-nums[i]]);
//                result.push_back(i);
//            }
//            saved.insert(make_pair(nums[i], i));
//        }
//        return result;
//    }
//};
//int main()
//{
//    map<int, int> m;
//    map<int, int>::iterator m_it;
//    m.insert(make_pair(1, 1));
//    m.insert(make_pair(4, 2));
//    m.insert(make_pair(2, 3));
//    for(m_it = m.begin(); m_it != m.end(); m_it++)
//    {
//        cout<< m_it->second;
//    }
//}
