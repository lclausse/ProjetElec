
int x1 = 90;
int y1 = 90;
PVector z1;

int x2 = 400;
int y2 = 80;
PVector z2;

int x3 = 350;
int y3 = 400;
PVector z3;

int xtarget = 300;
int ytarget = 300;
PVector ztarget;

int noise = 20;


void setup()
{
  size(1000, 1000);
  z1 = new PVector(x1, y1);
  z2 = new PVector(x2, y2);
  z3 = new PVector(x3, y3);
  ztarget = new PVector(xtarget, ytarget);
}


void draw()
{
  background(255);
  drawCircles(z1);
  drawCircles(z2);
  drawCircles(z3);

  text("P1", z1.x+5, z1.y);
  text("P2", z2.x+5, z2.y);
  text("P3", z3.x+5, z3.y);

  ztarget.x = mouseX;
  ztarget.y = mouseY;

  pushMatrix();
  stroke(0);
  strokeWeight(1);
  line(z1.x, z1.y, ztarget.x, ztarget.y);
  line(z2.x, z2.y, ztarget.x, ztarget.y);
  line(z3.x, z3.y, ztarget.x, ztarget.y);
  popMatrix();
}


void drawCircles(PVector z)
{
  pushMatrix();
  fill(0);
  noStroke();
  translate(z.x, z.y);
  ellipse(0, 0, 10, 10);

  int sizeToTarget = int(PVector.dist(z, ztarget));
  noFill();
  strokeWeight(noise);
  stroke(0, 100);
  ellipse(0, 0, 2*sizeToTarget, 2*sizeToTarget);
  popMatrix();
}
