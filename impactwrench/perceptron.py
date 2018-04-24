import numpy as np
import tensorflow as tf
import matplotlib.pyplot as plt


def perceptron(x,y):
    mu, sigma = x.mean(axis=0), x.std(axis=0)
    X_train = (X_train - mu) / sigma
    X_test = (X_test - mu) / sigma

    plt.scatter(X_train[:,0] ,y_train[:,0], label='class 0', marker='o')
    plt.scatter(X_train[:,1] ,y_train[:,0], label='class 1', marker='s')
    plt.xlabel('feature 1')
    plt.ylabel('feature 2')
    plt.legend()
    plt.show()

    g = tf.Graph()


    n_features = X_train.shape[1]

    with g.as_default() as g:

        # initialize model parameters
        features = tf.placeholder(dtype=tf.float32, 
                                  shape=[None, n_features], name='features')
        targets = tf.placeholder(dtype=tf.float32, 
                                 shape=[None, 1], name='targets')
        params = {
            'weights': tf.Variable(tf.zeros(shape=[n_features, 1], 
                                            dtype=tf.float32), name='weights'),
            'bias': tf.Variable([[0.]], dtype=tf.float32, name='bias')}

        # forward pass
        linear = tf.matmul(features, params['weights']) + params['bias']
        ones = tf.ones(shape=tf.shape(linear)) 
        zeros = tf.zeros(shape=tf.shape(linear))
        prediction = tf.where(tf.less(linear, 0.), zeros, ones, name='prediction')

        # weight update
        diff = targets - prediction
        weight_update = tf.assign_add(params['weights'], 
                                      tf.reshape(diff * features, (n_features, 1)))
        bias_update = tf.assign_add(params['bias'], diff)

        saver = tf.train.Saver()

    with tf.Session(graph=g) as sess:

        sess.run(tf.global_variables_initializer())

        i = 0
        for example, target in zip(x, y):
            feed_dict = {features: example.reshape(-1, n_features),
                         targets: target.reshape(-1, 1)}
            _, _ = sess.run([weight_update, bias_update], feed_dict=feed_dict)
            
            i += 1
            if i >= 4:
                break
            

        modelparams = sess.run(params)    
        print('Model parameters:\n', modelparams)

        saver.save(sess, "./my_perceptron.ckpt")

        pred = sess.run(prediction, feed_dict={features: x})
        errors = np.sum(pred.reshape(-1) != y_train)
        print('Number of training errors:', errors)

    import os
    with tf.Session(graph=g) as sess:
        saver.restore(sess, "./my_perceptron.ckpt" )

    for epoch in range(1):
        for example, target in zip(x, y):
            feed_dict = {features: example.reshape(-1, n_features),
                         targets: target.reshape(-1, 1)}
            _, _ = sess.run([weight_update, bias_update], feed_dict=feed_dict)
            modelparams = sess.run(params)

    saver.save(sess, "./my_perceptron.ckpt")

    pred = sess.run(prediction, feed_dict={features: x})
    train_errors = np.sum(pred.reshape(-1) != y_train)
    pred = sess.run(prediction, feed_dict={features: x})
    test_errors = np.sum(pred.reshape(-1) != y)

    print('Number of training errors', train_errors)
    print('Number of test errors', test_errors)

    x_min = -2
    y_min = ( -(modelparams['weights'][0] * x_min) / modelparams['weights'][1]
            -(modelparams['bias'] / modelparams['weights'][1]) )

    x_max = 2
    y_max = ( -(modelparams['weights'][0] * x_max) / modelparams['weights'][1]
            -(modelparams['bias'] / modelparams['weights'][1]) )


    fig, ax = plt.subplots(1, 2, sharex=True, figsize=(7, 3))

    ax[0].plot([x_max], [ y_max])
    ax[1].plot([x_min], [y_min])

    ax[0].scatter(x[:,0] ,y[:,0], label='class 0', marker='o')
    ax[0].scatter(x[:,1] ,y[:,0], label='class 1', marker='s')

    ax[1].scatter(x[:,0] ,y[:,0], label='class 0', marker='o')
    ax[1].scatter(x[:,0] ,y[:,0], label='class 0', marker='o')

    ax[1].legend(loc='upper left')
    plt.show()

