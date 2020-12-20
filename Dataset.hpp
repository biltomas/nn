#pragma once

#include <string>
#include <iostream>
#include <fstream>
#include <set>
#include "RowVector.hpp"

template <typename T>
struct item {
    RowVector<T> data;
    int label;
};

class DataLoader {
    const std::string vec_path;
    const std::string label_path;
    const std::string delimiter = ",";
    const size_t vector_size = 28 * 28;
    const std::set<char> possible_labels {'0', '1', '2', '3', '4', '5', '6', '7', '8', '9'};

    RowVector<float> parse_vector(const std::string& str_vec) const {
        size_t prev_position = 0;
        std::vector<float> data;
        for (size_t i = 0; i < vector_size; i++) {
            //std::cout << i;
            size_t current_position = str_vec.find(delimiter, prev_position); 
            //std::cout << "Start: " << prev_position << " End: " << current_position << std::endl;
            if (current_position == std::string::npos && i != vector_size - 1) {
                throw std::invalid_argument("Unexpected end of line");
            }
            const std::string str_num = str_vec.substr(prev_position, current_position - prev_position);
            prev_position = current_position + 1;
            //std::cout << "\"" << str_num << "\"" << std::endl;
            float num = std::stof(str_num);
            data.push_back(num);
        }
        if (data.size() != vector_size) {
            throw std::invalid_argument("Incorrect vector size");
        }
        return RowVector<float>(data);
    }

    int parse_label(const std::string& str_label) const {
        if (possible_labels.find(str_label[0]) == possible_labels.end() || str_label[1] == '\n') {
            throw std::invalid_argument("Invalid label");
        }
        return std::stoi(str_label);
    }

    item<float> parse_single(std::ifstream& vec_file, std::ifstream& label_file) const {
        char vec_buffer[4096];
        char label_buffer[3];
        vec_file.getline(vec_buffer, 4096);
        std::string raw_vec(vec_buffer);
        //std::cout << "Raw vec: " << raw_vec << std::endl;
        label_file.getline(label_buffer, 3);
        std::string raw_label(label_buffer);
        if (vec_file.fail()) {
            throw std::out_of_range("Possible buffer overflow for vec file. Increase the buffer size (1024)");
        }
        if (label_file.fail()) {
            throw std::out_of_range("Possible buffer overflow for label file. Increase the buffer size (1024)");
        }
        return item<float> {parse_vector(raw_vec), parse_label(raw_label)};
    }

public:
    DataLoader(const std::string vec_path, const std::string label_path)
        : vec_path(vec_path)
        , label_path(label_path) {};


    std::vector<item<float>> load() const {
        std::ifstream vec_file(vec_path);
        if (vec_file.fail() || !vec_file.is_open()) {
            throw std::invalid_argument("Vec file error");
        }
        std::ifstream label_file(label_path);
        if (label_file.fail() || !label_file.is_open()) {
            throw std::invalid_argument("Label file error");
        }
        std::vector<item<float>> data;
        while (!vec_file.eof() && !label_file.eof()) {
            try {
                data.push_back(parse_single(vec_file, label_file));
            } catch (const std::out_of_range& e) {
                if (vec_file.eof() && label_file.eof()) {
                    break;
                }
                throw e;
            }
        }
        return data;
    }
};
