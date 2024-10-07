export module bool3;

export class bool3
{
 public:
  bool b[3];
  struct
  {
    bool x, y, z;
  };

  bool3();
  bool3(bool x, bool y, bool z) : x(x), y(y), z(z) {};
};