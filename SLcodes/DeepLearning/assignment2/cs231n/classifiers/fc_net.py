import numpy as np

from cs231n.layers import *
from cs231n.layer_utils import *


class TwoLayerNet(object):
  """
  A two-layer fully-connected neural network with ReLU nonlinearity and
  softmax loss that uses a modular layer design. We assume an input dimension
  of D, a hidden dimension of H, and perform classification over C classes.
  
  The architecture should be affine - relu - affine - softmax.

  Note that this class does not implement gradient descent; instead, it
  will interact with a separate Solver object that is responsible for running
  optimization.

  The learnable parameters of the model are stored in the dictionary
  self.params that maps parameter names to numpy arrays.
  """
  
  def __init__(self, input_dim=3*32*32, hidden_dim=100, num_classes=10,
               weight_scale=1e-3, reg=0.0):
    """
    Initialize a new network.

    Inputs:
    - input_dim: An integer giving the size of the input
    - hidden_dim: An integer giving the size of the hidden layer
    - num_classes: An integer giving the number of classes to classify
    - dropout: Scalar between 0 and 1 giving dropout strength.
    - weight_scale: Scalar giving the standard deviation for random
      initialization of the weights.
    - reg: Scalar giving L2 regularization strength.
    """
    self.params = {}
    self.reg = reg
    
    ############################################################################
    # TODO: Initialize the weights and biases of the two-layer net. Weights    #
    # should be initialized from a Gaussian with standard deviation equal to   #
    # weight_scale, and biases should be initialized to zero. All weights and  #
    # biases should be stored in the dictionary self.params, with first layer  #
    # weights and biases using the keys 'W1' and 'b1' and second layer weights #
    # and biases using the keys 'W2' and 'b2'.                                 #
    ############################################################################
    self.params['W1'] = weight_scale * np.random.randn(input_dim, hidden_dim)
    self.params['b1'] = np.zeros(hidden_dim)
    self.params['W2'] = weight_scale * np.random.randn(hidden_dim, num_classes)
    self.params['b2'] = np.zeros(num_classes)
    ############################################################################
    #                             END OF YOUR CODE                             #
    ############################################################################


  def loss(self, X, y=None):
    """
    Compute loss and gradient for a minibatch of data.

    Inputs:
    - X: Array of input data of shape (N, d_1, ..., d_k)
    - y: Array of labels, of shape (N,). y[i] gives the label for X[i].

    Returns:
    If y is None, then run a test-time forward pass of the model and return:
    - scores: Array of shape (N, C) giving classification scores, where
      scores[i, c] is the classification score for X[i] and class c.

    If y is not None, then run a training-time forward and backward pass and
    return a tuple of:
    - loss: Scalar value giving the loss
    - grads: Dictionary with the same keys as self.params, mapping parameter
      names to gradients of the loss with respect to those parameters.
    """  
    scores = None
    ############################################################################
    # TODO: Implement the forward pass for the two-layer net, computing the    #
    # class scores for X and storing them in the scores variable.              #
    ############################################################################
    
    #extract parameters from dictionary
    W1 = self.params['W1']
    b1 = self.params['b1']
    W2 = self.params['W2']
    b2 = self.params['b2']
    reg = self.reg
    
    #forward pass
    x1, cache1 = affine_forward(X,W1,b1)  
    x2, cache2 = relu_forward(x1)
    scores, cache3 = affine_forward(x2,W2,b2)  

    ############################################################################
    #                             END OF YOUR CODE                             #
    ############################################################################

    # If y is None then we are in test mode so just return scores
    if y is None:
      return scores
    
    loss, grads = 0, {}
    ############################################################################
    # TODO: Implement the backward pass for the two-layer net. Store the loss  #
    # in the loss variable and gradients in the grads dictionary. Compute data #
    # loss using softmax, and make sure that grads[k] holds the gradients for  #
    # self.params[k]. Don't forget to add L2 regularization!                   #
    #                                                                          #
    # NOTE: To ensure that your implementation matches ours and you pass the   #
    # automated tests, make sure that your L2 regularization includes a factor #
    # of 0.5 to simplify the expression for the gradient.                      #
    ############################################################################
      
    #get scores from training labels y to calculate loss
    loss, dscores = softmax_loss(scores,y)

    #L2 regularization
    loss =  loss + 0.5*reg*(np.sum(W1*W1)+np.sum(W2*W2))

    #Calculate gradients
    #remind architecture : [affine (W1X) - relu(H=max(0,x))] - [affine(O=W2H) - 
    #                       softmax(1/1+e(-o))]    
    dx2,dw2,db2 = affine_backward(dscores,cache3)
    dx1 = relu_backward(dx2,cache2)
    dx,dw1,db1 = affine_backward(dx1,cache1)
    
    #gradients
    grads['W1'] = dw1 + reg*W1
    grads['b1'] = db1
    grads['W2'] = dw2 + reg*W2
    grads['b2'] = db2

    ############################################################################
    #                             END OF YOUR CODE                             #
    ############################################################################

    return loss, grads


class FullyConnectedNet(object):
  """
  A fully-connected neural network with an arbitrary number of hidden layers,
  ReLU nonlinearities, and a softmax loss function. This will also implement
  dropout and batch normalization as options. For a network with L layers,
  the architecture will be
  
  {affine - [batch norm] - relu - [dropout]} x (L - 1) - affine - softmax
  
  where batch normalization and dropout are optional, and the {...} block is
  repeated L - 1 times.
  
  Similar to the TwoLayerNet above, learnable parameters are stored in the
  self.params dictionary and will be learned using the Solver class.
  """

  def __init__(self, hidden_dims, input_dim=3*32*32, num_classes=10,
               dropout=0, use_batchnorm=False, reg=0.0,
               weight_scale=1e-2, dtype=np.float32, seed=None):
    """
    Initialize a new FullyConnectedNet.
    
    Inputs:
    - hidden_dims: A 
    of integers giving the size of each hidden layer.
    - input_dim: An integer giving the size of the input.
    - num_classes: An integer giving the number of classes to classify.
    - dropout: Scalar between 0 and 1 giving dropout strength. If dropout=0 then
      the network should not use dropout at all.
    - use_batchnorm: Whether or not the network should use batch normalization.
    - reg: Scalar giving L2 regularization strength.
    - weight_scale: Scalar giving the standard deviation for random
      initialization of the weights.
    - dtype: A numpy datatype object; all computations will be performed using
      this datatype. float32 is faster but less accurate, so you should use
      float64 for numeric gradient checking.
    - seed: If not None, then pass this random seed to the dropout layers. This
      will make the dropout layers deteriminstic so we can gradient check the
      model.
    """
    self.use_batchnorm = use_batchnorm #check if batch normalization
    self.use_dropout = dropout > 0 #check if dropout
    self.reg = reg
    self.num_layers = 1 + len(hidden_dims) # nb of hidden + output
    self.dtype = dtype
    self.params = {}    

    ############################################################################
    # TODO: Initialize the parameters of the network, storing all values in    #
    # the self.params dictionary. Store weights and biases for the first layer #
    # in W1 and b1; for the second layer use W2 and b2, etc. Weights should be #
    # initialized from a normal distribution with standard deviation equal to  #
    # weight_scale and biases should be initialized to zero.                   #
    #                                                                          #
    # When using batch normalization, store scale and shift parameters for the #
    # first layer in gamma1 and beta1; for the second layer use gamma2 and     #
    # beta2, etc. Scale parameters should be initialized to one and shift      #
    # parameters should be initialized to zero.                                #
    ############################################################################

    #initialize weights and biases for affine operation 
    #at hidden layer H1
    self.params['W1'] = weight_scale * np.random.randn(input_dim, hidden_dims[0])
    self.params['b1'] = np.zeros(hidden_dims[0])

    #initialize weights and biases affine operation 
    #at subsequent hidden layers      
    for i in range(0,self.num_layers-2):      

      #initialize W and b for affine operation      
      self.params["W" + str(i+2)] = weight_scale * np.random.randn(hidden_dims[i], hidden_dims[i+1]) #W2,W3 ....    
      self.params["b" + str(i+2)] = np.zeros(hidden_dims[i+1]) #b2,b3 ....

    #initialize gamma and beta when batch normalization 
    #between affine and ReLU operations
    #one pair for each layer unit
    if self.use_batchnorm:
        num_hidden = len(hidden_dims)
        for i in range(0,num_hidden):
        
          self.params["gamma" + str(i+1)] = np.ones(hidden_dims[i]) #gamma1, gamma2 ...         
          self.params["beta" + str(i+1)] = np.zeros(hidden_dims[i]) #beta1, beta2 ...         
           
    #initialize weights and biases for affine operation at output
    #no batch normalization after output layer
    self.params["W" + str(len(hidden_dims)+1)] = weight_scale * np.random.randn(hidden_dims[len(hidden_dims)-1],num_classes)      
    self.params["b" + str(len(hidden_dims)+1)] = np.zeros(num_classes)               
    
   
    ############################################################################
    #                             END OF YOUR CODE                             #
    ############################################################################

    # When using dropout we need to pass a dropout_param dictionary to each
    # dropout layer so that the layer knows the dropout probability and the mode
    # (train / test). You can pass the same dropout_param to each dropout layer.
    self.dropout_param = {}
    if self.use_dropout:
      self.dropout_param = {'mode': 'train', 'p': dropout}
      if seed is not None:
        self.dropout_param['seed'] = seed
    
    # With batch normalization we need to keep track of running means and
    # variances, so we need to pass a special bn_param object to each batch
    # normalization layer. You should pass self.bn_params[0] to the forward pass
    # of the first batch normalization layer, self.bn_params[1] to the forward
    # pass of the second batch normalization layer, etc.
    self.bn_params = []
    if self.use_batchnorm:
      self.bn_params = [{'mode': 'train'} for i in xrange(self.num_layers - 1)]
    
    # Cast all parameters to the correct datatype
    for k, v in self.params.iteritems():
      self.params[k] = v.astype(dtype)


  def loss(self, X, y=None):
    """
    Compute loss and gradient for the fully-connected net.

    Input / output: Same as TwoLayerNet above.
    """
    X = X.astype(self.dtype)
    mode = 'test' if y is None else 'train'

    # Set train/test mode for batchnorm params and dropout param since they
    # behave differently during training and testing.
    if self.dropout_param is not None:
      self.dropout_param['mode'] = mode   
    if self.use_batchnorm:
      for bn_param in self.bn_params:
        bn_param[mode] = mode

    scores = None
    ############################################################################
    # TODO: Implement the forward pass for the fully-connected net, computing  #
    # the class scores for X and storing them in the scores variable.          #
    #                                                                          #
    # When using dropout, you'll need to pass self.dropout_param to each       #
    # dropout forward pass.                                                    #
    #                                                                          #
    # When using batch normalization, you'll need to pass self.bn_params[0] to #
    # the forward pass for the first batch normalization layer, pass           #
    # self.bn_params[1] to the forward pass for the second batch normalization #
    # layer, etc.                                                              #
    ############################################################################
    
    #forward pass from input data through hidden layers         
    #init
    caches    = []
    rcaches   = []
    bn_caches = []
    dcaches   = []
    prev      =  X

    #loop through hidden layers
    for i in range(0,self.num_layers-1):             

      #affine operation
      #store its inputs for backward pass
      prev,cache = affine_forward(prev,self.params["W" + str(i+1)],self.params["b" + str(i+1)])                    
      caches.append(cache)      

      #when batch normalization
      #store its inputs for backward pass
      if self.use_batchnorm:          
          prev, cache = batchnorm_forward(prev, self.params["gamma" + str(i+1)], self.params["beta" + str(i+1)], self.bn_params[i])
          bn_caches.append(cache)
         
      #ReLU operation
      #store its inputs for backward pass
      prev,cache = relu_forward(prev)                        
      rcaches.append(cache)
      
      #when dropout operation
      if self.use_dropout:
          prev, cache = dropout_forward(prev,self.dropout_param)  
          dcaches.append(cache)            
            
    #forward pass output layer            
    scores, cacheOut = affine_forward(prev,self.params["W"+str(self.num_layers)],self.params["b" + str(self.num_layers)])

    # ############################################################################
    # #                             END OF YOUR CODE                             #
    # ############################################################################

    # If test mode return early
    if mode == 'test':
      return scores

    loss, grads = 0.0, {}
    # ############################################################################
    # # TODO: Implement the backward pass for the fully-connected net. Store the #
    # # loss in the loss variable and gradients in the grads dictionary. Compute #
    # # data loss using softmax, and make sure that grads[k] holds the gradients #
    # # for self.params[k]. Don't forget to add L2 regularization!               #
    # #                                                                          #
    # # When using batch normalization, you don't need to regularize the scale   #
    # # and shift parameters.                                                    #
    # #                                                                          #
    # # NOTE: To ensure that your implementation matches ours and you pass the   #
    # # automated tests, make sure that your L2 regularization includes a factor #
    # # of 0.5 to simplify the expression for the gradient.                      #
    # ############################################################################

    #softmax conversion of output scores to probabilities
    #loss evaluation, and dloss/dscores gradient at output
    loss, dscores = softmax_loss(scores,y)

    #Loss L2-regularization
    #sum weights^2 and regularize loss
    sumWeights = 0
    for i in range(0,self.num_layers):      
      Wi = self.params["W"+str(i+1)]
      sumWeights += np.sum(Wi*Wi)                
    loss =  loss + 0.5*self.reg*sumWeights    

    #Backward pass gradients at output
    dx3,dw3,db3 = affine_backward(dscores,cacheOut)      
    grads["W"+str(self.num_layers)] = dw3 + self.reg*self.params["W"+str(self.num_layers)]
    grads["b"+str(self.num_layers)] = db3
    
    #gradients at each hidden layer           
    dup = dx3
    for i in range(0,self.num_layers-1):      
      #go backward 
      j   = self.num_layers-(i+2)    
      
      #when dropout operation
      if self.use_dropout:
          dup = dropout_backward(dup,dcaches[j]) 
    
      #ReLU  
      dup = relu_backward(dup,rcaches[j])    
        
      #when batch normalization with its gradients (w/o regularization)
      if self.use_batchnorm:          
          dup, dgammaj, dbetaj = batchnorm_backward(dup,bn_caches[j])
          grads["gamma" + str(j+1)] = dgammaj
          grads["beta" + str(j+1)] = dbetaj        
        
      #affine and its gradients (w/o regularization)
      dup,dwj,dbj = affine_backward(dup,caches[j])          
      grads["W" + str(j+1)] = dwj + self.reg*self.params["W" + str(j+1)]
      grads["b" + str(j+1)] = dbj           
             

    # ############################################################################
    # #                             END OF YOUR CODE                             #
    # ############################################################################

    return loss, grads
