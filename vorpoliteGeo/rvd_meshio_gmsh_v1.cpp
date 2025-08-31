#include "rvd_meshio_gmsh_v1.h"

#include <fstream>
#include <sstream>
#include <unordered_map>
#include <vector>
#include <string>
#include <cctype>

namespace rvd {

    // Internal utility functions - completely isolated
    namespace loader_utils {

        std::string remove_whitespace(const std::string& input) {
            std::size_t start = 0;
            std::size_t end = input.size();

            while (start < end && std::isspace(static_cast<unsigned char>(input[start]))) {
                ++start;
            }

            while (end > start && std::isspace(static_cast<unsigned char>(input[end - 1]))) {
                --end;
            }

            return input.substr(start, end - start);
        }

        bool strings_equal_ignore_case(const std::string& a, const char* b) {
            std::size_t b_len = 0;
            while (b[b_len] != '\0') ++b_len; // Manual strlen to avoid macro issues

            if (a.size() != b_len) return false;

            for (std::size_t i = 0; i < a.size(); ++i) {
                char a_char = static_cast<char>(std::tolower(static_cast<unsigned char>(a[i])));
                char b_char = static_cast<char>(std::tolower(static_cast<unsigned char>(b[i])));
                if (a_char != b_char) return false;
            }
            return true;
        }

        int get_element_node_count(int element_type) {
            switch (element_type) {
            case 2: return 3; // 3-node triangle
            case 4: return 4; // 4-node tetrahedron
            default: return -1; // unsupported
            }
        }
    }

    bool load_gmsh_v1_ascii(const std::string& path, TetMesh& output, std::string* error_message) {
        std::ifstream file_stream(path.c_str());
        if (!file_stream.is_open()) {
            if (error_message) *error_message = "Cannot open file: " + path;
            return false;
        }

        output.points.clear();
        output.tets.clear();
        output.tris.clear();

        std::unordered_map<long long, std::uint32_t> node_id_to_index;

        std::string current_line;
        while (std::getline(file_stream, current_line)) {
            current_line = loader_utils::remove_whitespace(current_line);

            if (current_line.empty()) continue;
            if (current_line[0] == '#') continue;
            if (current_line.size() >= 2 && current_line[0] == '/' && current_line[1] == '/') continue;

            // Handle nodes section
            if (loader_utils::strings_equal_ignore_case(current_line, "$NOD") ||
                loader_utils::strings_equal_ignore_case(current_line, "$Nodes")) {

                bool is_legacy_format = loader_utils::strings_equal_ignore_case(current_line, "$NOD");

                if (!std::getline(file_stream, current_line)) {
                    if (error_message) *error_message = "Unexpected EOF after $NOD/$Nodes";
                    return false;
                }

                std::istringstream line_parser(loader_utils::remove_whitespace(current_line));
                std::size_t node_count = 0;
                if (!(line_parser >> node_count)) {
                    if (error_message) *error_message = "Could not read node count";
                    return false;
                }

                output.points.clear();
                output.points.resize(node_count);
                node_id_to_index.clear();
                node_id_to_index.reserve(node_count);

                for (std::size_t node_index = 0; node_index < node_count; ++node_index) {
                    if (!std::getline(file_stream, current_line)) {
                        if (error_message) *error_message = "Unexpected EOF in nodes block";
                        return false;
                    }

                    current_line = loader_utils::remove_whitespace(current_line);
                    if (current_line.empty()) {
                        --node_index;
                        continue;
                    }

                    std::istringstream node_parser(current_line);
                    long long node_id = 0;
                    double x_coord = 0.0, y_coord = 0.0, z_coord = 0.0;

                    if (!(node_parser >> node_id >> x_coord >> y_coord >> z_coord)) {
                        if (error_message) *error_message = "Bad node line: " + current_line;
                        return false;
                    }

                    std::uint32_t array_index = static_cast<std::uint32_t>(node_index);
                    node_id_to_index[node_id] = array_index;

                    output.points[array_index][0] = x_coord;
                    output.points[array_index][1] = y_coord;
                    output.points[array_index][2] = z_coord;
                }

                if (!std::getline(file_stream, current_line)) {
                    if (error_message) *error_message = "Missing end of nodes block";
                    return false;
                }
                current_line = loader_utils::remove_whitespace(current_line);

                bool valid_end_token = false;
                if (is_legacy_format) {
                    valid_end_token = loader_utils::strings_equal_ignore_case(current_line, "$ENDNOD") ||
                        loader_utils::strings_equal_ignore_case(current_line, "$EndNodes");
                }
                else {
                    valid_end_token = loader_utils::strings_equal_ignore_case(current_line, "$EndNodes") ||
                        loader_utils::strings_equal_ignore_case(current_line, "$ENDNOD");
                }

                if (!valid_end_token) {
                    if (error_message) *error_message = "Expected $ENDNOD/$EndNodes";
                    return false;
                }
                continue;
            }

            // Handle elements section
            if (loader_utils::strings_equal_ignore_case(current_line, "$ELM") ||
                loader_utils::strings_equal_ignore_case(current_line, "$Elements")) {

                bool is_legacy_format = loader_utils::strings_equal_ignore_case(current_line, "$ELM");

                if (!std::getline(file_stream, current_line)) {
                    if (error_message) *error_message = "Unexpected EOF after $ELM/$Elements";
                    return false;
                }

                std::istringstream line_parser(loader_utils::remove_whitespace(current_line));
                std::size_t element_count = 0;
                if (!(line_parser >> element_count)) {
                    if (error_message) *error_message = "Could not read element count";
                    return false;
                }

                output.tets.clear();
                output.tris.clear();
                output.tets.reserve(element_count);
                output.tris.reserve(element_count);

                for (std::size_t element_index = 0; element_index < element_count; ++element_index) {
                    if (!std::getline(file_stream, current_line)) {
                        if (error_message) *error_message = "Unexpected EOF in elements block";
                        return false;
                    }

                    current_line = loader_utils::remove_whitespace(current_line);
                    if (current_line.empty()) {
                        --element_index;
                        continue;
                    }

                    std::istringstream element_parser(current_line);
                    std::vector<long long> tokens;
                    tokens.reserve(16);

                    long long token_value;
                    while (element_parser >> token_value) {
                        tokens.push_back(token_value);
                    }

                    if (tokens.size() < 3) {
                        if (error_message) *error_message = "Bad element line: " + current_line;
                        return false;
                    }

                    int element_type = static_cast<int>(tokens[1]);
                    int required_nodes = loader_utils::get_element_node_count(element_type);

                    if (required_nodes < 0) continue; // skip unsupported element types

                    if (static_cast<int>(tokens.size()) < 2 + required_nodes) {
                        if (error_message) *error_message = "Element line too short for node list: " + current_line;
                        return false;
                    }

                    // Extract node IDs from the end of the token list
                    int node_start_index = static_cast<int>(tokens.size()) - required_nodes;

                    if (element_type == 4) { // Tetrahedron
                        std::array<std::uint32_t, 4> tet_nodes{};

                        for (int i = 0; i < 4; ++i) {
                            long long node_id = tokens[node_start_index + i];
                            auto lookup = node_id_to_index.find(node_id);
                            if (lookup == node_id_to_index.end()) {
                                if (error_message) {
                                    *error_message = "Unknown node id in tet: ";
                                    // Manual number to string conversion to avoid macro issues
                                    std::ostringstream converter;
                                    converter << node_id;
                                    *error_message += converter.str();
                                }
                                return false;
                            }
                            tet_nodes[i] = lookup->second;
                        }
                        output.tets.push_back(tet_nodes);
                    }
                    else if (element_type == 2) { // Triangle
                        std::array<std::uint32_t, 3> tri_nodes{};

                        for (int i = 0; i < 3; ++i) {
                            long long node_id = tokens[node_start_index + i];
                            auto lookup = node_id_to_index.find(node_id);
                            if (lookup == node_id_to_index.end()) {
                                if (error_message) {
                                    *error_message = "Unknown node id in tri: ";
                                    // Manual number to string conversion to avoid macro issues
                                    std::ostringstream converter;
                                    converter << node_id;
                                    *error_message += converter.str();
                                }
                                return false;
                            }
                            tri_nodes[i] = lookup->second;
                        }
                        output.tris.push_back(tri_nodes);
                    }
                }

                if (!std::getline(file_stream, current_line)) {
                    if (error_message) *error_message = "Missing end of elements block";
                    return false;
                }
                current_line = loader_utils::remove_whitespace(current_line);

                bool valid_end_token = false;
                if (is_legacy_format) {
                    valid_end_token = loader_utils::strings_equal_ignore_case(current_line, "$ENDELM") ||
                        loader_utils::strings_equal_ignore_case(current_line, "$EndElements");
                }
                else {
                    valid_end_token = loader_utils::strings_equal_ignore_case(current_line, "$EndElements") ||
                        loader_utils::strings_equal_ignore_case(current_line, "$ENDELM");
                }

                if (!valid_end_token) {
                    if (error_message) *error_message = "Expected $ENDELM/$EndElements";
                    return false;
                }
                continue;
            }

            // Ignore all other sections
        }

        if (output.points.empty() || output.tets.empty()) {
            if (error_message) *error_message = "File parsed, but points or tets are empty";
            return false;
        }

        return true;
    }

} // namespace rvd