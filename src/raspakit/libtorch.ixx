module;

#include <iostream>
#include <iomanip>
#ifdef BUILD_LIBTORCH
#include <torch/torch.h>
#endif

export module libtorch_test;

#ifdef BUILD_LIBTORCH

export void test_libtorch() 
{
   std::cout << "Linear Regression\n\n";

    // Device
    auto mps_available = torch::mps::is_available();
    torch::Device device(mps_available ? torch::kMPS : torch::kCPU);
    std::cout << (mps_available ? "Metal Performance Shaders available. Training on GPU." : "Training on CPU.") << '\n';

    // Hyper parameters
    const int64_t input_size = 1;
    const int64_t output_size = 1;
    const size_t num_epochs = 60;
    const double learning_rate = 0.001;

    // Sample dataset
    auto x_train = torch::randint(0, 10, {15, 1},
                                  torch::TensorOptions(torch::kFloat).device(device));

    auto y_train = torch::randint(0, 10, {15, 1},
                                  torch::TensorOptions(torch::kFloat).device(device));

    // Linear regression model
    torch::nn::Linear model(input_size, output_size);
    model->to(device);

    // Optimizer
    torch::optim::SGD optimizer(model->parameters(), torch::optim::SGDOptions(learning_rate));

    // Set floating point output precision
    std::cout << std::fixed << std::setprecision(4);

    std::cout << "Training...\n";

    // Train the model
    for (size_t epoch = 0; epoch != num_epochs; ++epoch) {
        // Forward pass
        auto output = model->forward(x_train);
        auto loss = torch::nn::functional::mse_loss(output, y_train);

        // Backward pass and optimize
        optimizer.zero_grad();
        loss.backward();
        optimizer.step();

        if ((epoch + 1) % 5 == 0) {
            std::cout << "Epoch [" << (epoch + 1) << "/" << num_epochs <<
                "], Loss: " << loss.item<double>() << "\n";
        }
    }

    std::cout << "Training finished!\n";
}
#endif
