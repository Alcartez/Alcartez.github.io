# Flask Tutorial

### Basic Flask App

```python
# Flask print basic variable.

from flask import Flask # Import Flask

app = Flask(__name__) # Assigning app variable

@app.route("/") # Routing it to the subpage in a domain , in this case it'll be the index page cuz no page name is defined.

def home(): # Function home to print the statement in HTML
    return "Hello! This is the main page <h1> Hello </h1>." 

@app.route("/<name>") # We assign a variable as a subpage name and when we execute it , the variable is stored in ram.
def user(name): # This stored variable name is used in this function.
    return f"Hello {name}"

#------------------------------------------#    
#  Something you should know is that Python assigns the name "__main__" to the script when the script is executed.
#  If the script is imported from another script, the script keeps it given name (e.g. hello.py).
#  In our case, we are executing the script. Therefore, __name__ will be equal to "__main__".
#  That means the if conditional statement is satisfied and the app.run() method will be executed. 
#  This technique allows the programmer to have control over the script’s behavior.
#------------------------------------------#    


if __name__ == "__main__": # Run the app using condition.
    app.run()  
```



### Redirect

Lets demonstrate the redirect function. Used to redirect user to another
page upon function execution.

```python
return redirect(url_for("user" , name = "ADMIN"))
```



```python
# REDIRECT #

from flask import Flask , redirect, url_for # Import Flask , redirect and url_for

app = Flask(__name__) # Assigning app variable

@app.route("/") # Routing it to the subpage in a domain , in this case it'll be the index page cuz no page name is defined.

def home(): # Function home to print the statement in HTML
    return "Hello! This is the main page <h1> Hello </h1>." 

@app.route("/<name>") # We assign a variable as a subpage name and when we execute it , the variable is stored in ram.
def user(name): # This stored variable name is used in this function.
    return f"Hello {name}"

@app.route("/admin") # If you use "/admin/" , you can access it with or without the end "/"
def admin():
    # return redirect(url_for("home"))   # Instead of redirecting to the page with "/" , we redirect it to the function i.e. home().
    return redirect(url_for("user" , name = "Admin!" )) # For redirecting to functions with variable inputs.

if __name__ == "__main__": # Run the app using condition.
    app.run()  
```



## HTML TEMPLATES

#### flask.py 

```python
def home(name):
    return render_template("index.html" , content=name) 
```

The content variable will be used in the html code below.

#### index.html 

```html
<!DOCTYPE html>
<html>
    <head>
        <title>Home Page</title>
    </head>
    <body>
        <h1>Home Page</h1>
        <p>{{content}}</p>
    </body>
</html>
```



```python
# HTML TEMPLATES #

# Start by creating a specific directory called "templates" inside the base directory.
# Create an html file of any name that you want to be the template page.
# Create {{content}} in the html file

from flask import Flask , redirect, url_for , render_template # Import Flask , redirect and url_for

app = Flask(__name__) # Assigning app variable

@app.route("/<name>") # Routing it to the subpage in a domain , in this case it'll be the index page cuz no page name is defined.
def home(name):
    return render_template("index.html" , content=name) # Content is a variable called in the html file
    # You can add as many variables as you want , as long as they are called. 
    
if __name__ == "__main__": # Run the app using condition.
    app.run()  
```



### You can write python code in the html files

#### example_loop_static.html

```html
<!DOCTYPE html>
<html>
    <head>
        <title>Home Page</title>
    </head>
    <body>
        <h1>Home Page</h1>
        {% for _ in range(10)%}
            <p>Hello</p>
        {% endfor %}       
    </body>
</html>
```

This will print hello 10 times

#### example_loop_variable.html 

```html
<!DOCTYPE html>
<html>
    <head>
        <title>Home Page</title>
    </head>
    <body>
        <h1>Home Page</h1>
        {% for x in range(10)%}
            <p>{{x}}</p>
        {% endfor %}       
    </body>
</html>
```

This will print the values 1 to 10

#### flask.py 

```python
def home(name):
    return render_template("index.html" , list=["1","2","3","4"]) 
```

#### example_loop_list.html 

```html
<!DOCTYPE html>
<html>
    <head>
        <title>Home Page</title>
    </head>
    <body>
        <h1>Home Page</h1>
        {% for x in list%}
            <p>{{x}}</p>
        {% endfor %}       
    </body>
</html>
```

This will print the values the values in the list

#### example_condition.html 

```html
<!DOCTYPE html>
<html>
    <head>
        <title>Home Page</title>
    </head>
    <body>
        <h1>Home Page</h1>
        {% raw %}
        {% for x in list%}
            {% if x == "3" %}
                <h1>THREE</h1>
            {% elif x == "2" %}
                <h2>TWO</h2>
            {% else %}
                <h3>{{x}}</h3>
            {% endif %}
        {% endfor %} 
        {% endraw %}
    </body>
</html>
```



## TEMPLATE INHERITANCE

### Parent

The parent html is used as a template for the child templates to give
out an uniform sense of page.

#### parent.html 

```html
<!DOCTYPE html>
<html>
    <head>
    <title>
        {% raw %}
        {% block content %}{% endblock %}
        {% endraw %}
    </title>
    </head>
    <body>
        <h1>WEBSITE</h1>
    </body>
</html>
```

### Child

The child html inherits the template from the parent html and uses the
same layout for displaying everything. The content blocks are used to
put content in the same positions as in the parent block.

#### child.html 

```html
{% raw %}
{% extends "parent.html" %}
{% block content %}Home Page{% endblock %}
{% block content %}
<h1> Test </h1>
{% endblock %}
{% endraw %}
```



```python
# TEMPLATE INHERITANCE #
from flask import Flask , redirect, url_for , render_template # Import Flask , redirect and url_for

app = Flask(__name__) 

@app.route("/parent") 
def parent():
    return render_template("parent.html")

@app.route("/child")
def child():
    return render_template("child.html") 

@app.route("/<page>/") 
def redirect(page):
    return redirect(url_for("parent")) 
    
    
if __name__ == "__main__": 
    app.run(debug=True)  # Activate debugger while running
```



## HTTP Methods (GET/POST)

#### login.html 

```html
{% extends "base.html" %}
{% block title %}Login Page{% endblock %}

{% block content %}
<form action = "#" method="post">
    <p>Name:<p>
    <p><input type="text" name="nm"></p>
    <p><input type="submit" value = "SUBMIT">
</form>
{% endblock %}
```



```python
from flask import Flask , redirect , url_for , render_template , request

app = Flask(__name__)

@app.route("/")
def home():
    return render_template("parent.html")

@app.route("/login", methods=["POST","GET"]) #Using POST and GET in this page.
def login():
    if request.method == "POST":
        user = request.form["nm"] # Request the name from the form.
        return redirect(url_for("user",usr=user)) # Redirect to other page
    else :
        return render_template("login.html")

@app.route("/<usr>")
def user(usr):
    return f"<h1>{usr}</h1>"

if __name__ == "__main__":
    app.run(debug=True)
```



# Sessions

#### login.html 

```html
{% extends "base.html" %}
{% block title %}Login Page{% endblock %}

{% block content %}
<form action = "#" method="post">
    <p>Name:<p>
    <p><input type="text" name="nm"></p>
    <p><input type="submit" value = "SUBMIT">
</form>
{% endblock %}
```



```python
# SESSIONS 

from flask import Flask , redirect , url_for , render_template , request , session
from datetime import timedelta

app = Flask(__name__)
app.secret_key = "hello" # For session encryption
app.permanent_session_lifetime = timedelta(days=5) # Store permanent session data for 5 days


@app.route("/")
def home():
    return render_template("parent.html")

@app.route("/login", methods=["POST","GET"]) #Using POST and GET in this page.
def login():
    if request.method == "POST":
        # session.permanent = True  # To make session permanent
        user = request.form["nm"] # Request the name from the form.
        session["user"] = user # Saving user in session
        return redirect(url_for("user",usr=user)) # Redirect to other page
    else :
        if "user" in session :
             return redirect(url_for("user"))
        
        return render_template("login.html")

@app.route("/user")
def user(usr):
    if "user" in session:
        user = session["user"] # Getting user from session
        return f"<h1>{usr}</h1>"
    else :
        return redirect(url_for("login")) # If nothing in session , redirect to login.

@app.route("/logout")
def logout():
    session.pop("user",none) # Close the session
    return redirect(url_for("login"))

if __name__ == "__main__":
    app.run(debug=True)
```



## Message Flashing

#### login.html 

```html
{% extends "base.html" %}
{% block title %}Login Page{% endblock %}

{% block content %}
    {% with messages = get_flashed_messages() %}
        {% if messages %}
            {% for msg in messages %}
                <p>{{msg}}</p>
            {% endfor %}
        {% endif %}
    {% endwith %}
    
<form action = "#" method="post">
    <p>Name:<p>
    <p><input type="text" name="nm"></p>
    <p><input type="submit" value = "SUBMIT">
</form>
{% endblock %}
```

#### user.html 

```html
{% extends "base.html" %}
{% block title %}User Page{% endblock %}

{% block content %}
    {% with messages = get_flashed_messages() %}
        {% if messages %}
            {% for msg in messages %}
                <p>{{msg}}</p>
            {% endfor %}
        {% endif %}
    {% endwith %}
<h1>Welcome {{user}}</h1>
{% endblock %}
```



```python
# Message Flashing

from flask import Flask , redirect , url_for , render_template , request , session , flash
from datetime import timedelta

app = Flask(__name__)
app.secret_key = "hello" # For session encryption
app.permanent_session_lifetime = timedelta(days=5) # Store permanent session data for 5 days


@app.route("/")
def home():
    return render_template("index.html")

@app.route("/login", methods=["POST","GET"]) #Using POST and GET in this page.
def login():
    if request.method == "POST":
        # session.permanent = True  # To make session permanent
        user = request.form["nm"] # Request the name from the form.
        session["user"] = user # Saving user in 
        flash(f"Successfully logged in ! {user}", "info")
        return redirect(url_for("user",usr=user)) # Redirect to other page
    else :
        if "user" in session :
            flash(f"Already logged in!")
            return redirect(url_for("user"))
        return render_template("login.html")

@app.route("/user")
def user(usr):
    if "user" in session:
        user = session["user"] # Getting user from session
        return render_template("user.html" , user=user)
    else :
        flash(f"Not logged in!")
        return redirect(url_for("login")) # If nothing in session , redirect to login.

@app.route("/logout")
def logout():
    session.pop("user",none) # Close the session
    if "user" in session:
        flash(f"You have been logged out ! {user}", "info")
    return redirect(url_for("login"))

if __name__ == "__main__":
    app.run(debug=True)
```



## SQLAlchemy

```python
# To install SqlAlchemy
pip install flask-sqlalchemy
```

#### parent.html 

```html
<!DOCTYPE html>
<html>
    <head>
    <title>{% block content %}{% endblock %}</title>
    </head>
    <body>
        <h1>WEBSITE</h1>
        <div class ="container-fluid">
            {% block content %}
            {% endblock %}
        </div>
    </body>
</html>
```

#### user.html 

```html
{% extends "parent.html" %}
{% block title %}User Page{% endblock %}

{% block content %}
    {% with messages = get_flashed_messages() %}
        {% if messages %}
            {% for msg in messages %}
                <p>{{msg}}</p>
            {% endfor %}
        {% endif %}
    {% endwith %}
<form action ='#' method="POST">
    <input type ="email" name = "email" placeholder = "Enter Email" value="{{email if email}}"/>
</form>
{% endblock %}
```



```python
# SQLAlchemy

from flask import Flask , redirect , url_for , render_template , request , session , flash
from datetime import timedelta
from flask_sqlalchemy import SQLAlchemy

app = Flask(__name__)
app.secret_key = "hello" # For session encryption
app.config['SQLALCHEMY_DATABASE_URI'] = 'sqlite:///users.sqlite3'
app.config['SQLALCHEMY_TRACK_MODIFICATIONS'] = False
app.permanent_session_lifetime = timedelta(minutes=5) # Store permanent session data for 5 days

db = SQLAlchemy(app)

class users(db.Model): # Make a database class for object
    _id = db.Column("id" ,db.Integer, primary_key=True) # Primary unique key
    name = db.Column(db.String(100))
    email = db.Column(db.String(100))

    def __init__(self, name , email):
        self.name = name
        self.email = email

@app.route("/")
def home():
    return render_template("index.html")

@app.route("/login", methods=["POST","GET"]) #Using POST and GET in this page.
def login():
    if request.method == "POST":
        # session.permanent = True  # To make session permanent
        user = request.form["nm"] # Request the name from the form.
        session["user"] = user # Saving user in 
        flash(f"Successfully logged in ! {user}", "info")
        return redirect(url_for("user",usr=user)) # Redirect to other page
    else :
        if "user" in session :
            flash(f"Already logged in!")
            return redirect(url_for("user"))
        return render_template("login.html")

@app.route("/user" , methods=["POST","GET"])
def user():
    email = None
    if "user" in session:
        user = session["user"] # Getting user from session
        if request.method == "POST":
            email = request.form["email"]
            session["email"] = email
            flash("Email was saved!")
        else:
            if "email" in session :
                email = session["email"]

        return render_template("user.html" , email=email)

    else :
        flash(f"Not logged in!")
        return redirect(url_for("login")) # If nothing in session , redirect to login.

@app.route("/logout")
def logout():
    session.pop("user",none) # Close the session
    if "user" in session:
        flash(f"You have been logged out ! {user}", "info")
    return redirect(url_for("login"))

if __name__ == "__main__":
    db.create_all() # Create Database
    app.run(debug=True)
```



# Adding, Deleting & Updating Users 
## SQLAlchemy 
```python
# To install SqlAlchemy
pip install flask-sqlalchemy
```

#### parent.html 

```html
<!DOCTYPE html>
<html>
    <head>
    <title>{% block content %}{% endblock %}</title>
    </head>
    <body>
        <h1>WEBSITE</h1>
        <div class ="container-fluid">
            {% block content %}
            {% endblock %}
        </div>
    </body>
</html>
```

#### user.html 

```html
{% extends "parent.html" %}
{% block title %}User Page{% endblock %}

{% block content %}
    {% with messages = get_flashed_messages() %}
        {% if messages %}
            {% for msg in messages %}
                <p>{{msg}}</p>
            {% endfor %}
        {% endif %}
    {% endwith %}
<form action ='#' method="POST">
    <input type ="email" name = "email" placeholder = "Enter Email" value="{{email if email}}"/>
</form>
{% endblock %}
```

#### view.html 

```html
{% extends "parent.html" %}
{% block title %} View All Users {% endblock %}
{% block content %}
    {% for item in values %}
        <p>Name : {{item.name}} , Email : {{item.email}}</p>
    {% endfor %}
{% endblock %}
```



```python
# SQLAlchemy Advanced

from flask import Flask , redirect , url_for , render_template , request , session , flash
from datetime import timedelta
from flask_sqlalchemy import SQLAlchemy

app = Flask(__name__)
app.secret_key = "hello" # For session encryption
app.config['SQLALCHEMY_DATABASE_URI'] = 'sqlite:///users.sqlite3'
app.config['SQLALCHEMY_TRACK_MODIFICATIONS'] = False
app.permanent_session_lifetime = timedelta(minutes=5) # Store permanent session data for 5 days

db = SQLAlchemy(app)

class users(db.Model): # Make a database class for object
    _id = db.Column("id" ,db.Integer, primary_key=True) # Primary unique key
    name = db.Column(db.String(100))
    email = db.Column(db.String(100))

    def __init__(self, name , email):
        self.name = name
        self.email = email

@app.route("/")
def home():
    return render_template("index.html")

@app.route("/view")
def view():
    return render_template("view.html" , values = users.query.add())

@app.route("/login", methods=["POST","GET"]) #Using POST and GET in this page.
def login():
    if request.method == "POST":
        # session.permanent = True  # To make session permanent
        user = request.form["nm"] # Request the name from the form.
        session["user"] = user # Saving user in 
        
        #delete_user = users.query.filter_by(name = user).delete() # Will be deleted when it's commited 
        #for user in found_user: # To delete multiple users
        #   user.delete()
             
        found_user = users.query.filter_by(name = user).first()
        if found_user:
            session["email"] = found_user.email

        else:
            usr = users(user , None) # Can replace None with ""
            db.session.add(usr) # Add the user
            db.session.commit() # Commit changes

        flash(f"Successfully logged in ! {user}", "info")
        return redirect(url_for("user",usr=user)) # Redirect to other page
    else :
        if "user" in session :
            flash(f"Already logged in!")
            return redirect(url_for("user"))
        return render_template("login.html")

@app.route("/user" , methods=["POST","GET"])
def user():
    email = None
    if "user" in session:
        user = session["user"] # Getting user from session
        if request.method == "POST":
            email = request.form["email"]
            session["email"] = email
            found_user = users.query.filter_by(name = user).first()
            found_user.email = email
            db.session.commit()
            flash("Email was saved!")
        else:
            if "email" in session :
                email = session["email"]

        return render_template("user.html" , email=email)

    else :
        flash(f"Not logged in!")
        return redirect(url_for("login")) # If nothing in session , redirect to login.

@app.route("/logout")
def logout():
    session.pop("user",none) # Close the session
    if "user" in session:
        flash(f"You have been logged out ! {user}", "info")
    return redirect(url_for("login"))

if __name__ == "__main__":
    db.create_all() # Create Database
    app.run(debug=True)
```



# Static Files

Note : Bootstrap link is to be above the css link.

#### base.html 
```html
<html>
    <head>
        <title>{% block title %} {% endblock %}</title>
        <link rel="stylesheet" type = "text/css" href="{{url_for('static', filename='style.css')}}">
    </head>
    <body>
        {% block content %}
        {% endblock %}
    </body>
</html>
```

#### home.html 

```html
{% extends 'base.html'%}
{% block title %}
Home Page
{% endblock %}
{% block content %}
<h1> Home </h1>
<image src="{{url_for('static', filename='image.jpeg')}}">
{% endblock %}
```

To be saved in \"base_dir/static/styles/style.css\"

#### style.css

``` {.css}
body {color:#fffa;}
```



```python
from flask import Flask , render_template

app = Flask(__name__)

@app.route("/home")
@app.route("/")
def home():
    return render_template("home.html")

if __name__ == "__main__":
    app.run(debug=True)
```



# Blue Prints



```python
# main.py

from flask import Flask , render_template
from second import second

app = Flask(__name__)
app.register_blueprint(second , url_prefix="")

@app.route("/")
def test():
    return "<h1>Test</h1>"

if __name__ == "__main__":
    app.run(debug=True)
```



```python
# second.py

from flask import Flask , render_template

second = Blueprint("second" , __name__ , static_folder="static" , template_folder = "templates")

@second.route("/home")
@second.route("/")
def home():
    return render_template("home.html")

@second.route("/test")
def test():
    return "<h1>Test_2b</h1>"
```

