use orbclient::{Color, EventOption, Renderer, Window, WindowFlag};
use std::{cmp, f64};

use self::vec::Vec2;
mod vec;

struct Circle {
    p: Vec2,
    r: f64,
}

impl Circle {
    fn new(p: Vec2, r: f64) -> Self {
        Self { p, r }
    }

    fn draw<R: Renderer>(&self, renderer: &mut R, color: Color, scale: f64) {
        let w = renderer.width() as i32;
        let h = renderer.height() as i32;
        renderer.circle(
            (self.p.x * scale).round() as i32 + w / 2,
            (self.p.y * scale).round() as i32 + h / 2,
            -cmp::max(1, (self.r * scale * 10.0).round() as i32),
            color
        );
    }

    fn intersection(&self, other: &Self) -> Option<(Vec2, f64)> {
        let d = other.p - self.p;
        let r = d.length();

        if r < (self.r + other.r) {
            Some((
                d,
                r
            ))
        } else {
            None
        }
    }
}

struct Line {
    a: Vec2,
    b: Vec2
}

impl Line {
    fn new(a: Vec2, b: Vec2) -> Self {
        Self { a, b }
    }

    fn draw<R: Renderer>(&self, renderer: &mut R, color: Color, scale: f64) {
        let w = renderer.width() as i32;
        let h = renderer.height() as i32;
        renderer.wu_line(
            (self.a.x * scale).round() as i32 + w / 2,
            (self.a.y * scale).round() as i32 + h / 2,
            (self.b.x * scale).round() as i32 + w / 2,
            (self.b.y * scale).round() as i32 + h / 2,
            color
        );
    }
}

struct Object {
    p: Vec2,
    v: Vec2,
    r: f64,
    m: f64,
    c: Color,
}

impl Object {
    fn collision_geometry(&self, t: f64) -> (Circle, Line, Line, Circle) {
        let n = self.p + self.v * t;
        let theta = self.v.angle();
        let current = Circle::new(self.p, self.r);
        let left = {
            let theta = theta - f64::consts::PI/2.0;
            let (sin, cos) = theta.sin_cos();
            let offset = Vec2::new(self.r * cos, self.r * sin);
            Line::new(self.p + offset, n + offset)
        };
        let right = {
            let theta = theta + f64::consts::PI/2.0;
            let (sin, cos) = theta.sin_cos();
            let offset = Vec2::new(self.r * cos, self.r * sin);
            Line::new(self.p + offset, n + offset)
        };
        let next = Circle::new(n, self.r);
        (current, left, right, next)
    }

    // Check for collision point
    fn collide(&self, other: &Self, t: f64) -> Option<(Vec2, Vec2)> {
        // When colliding objects, the following mut not collide:
        // - The current position of the object and it's circular radius
        // - The next position of the object and it's circular radius
        // - The rectangle from the outside of the radius extending in the direction of velocity

        let (self_current, self_left, self_right, self_next) = self.collision_geometry(t);
        let (other_current, other_left, other_right, other_next) = other.collision_geometry(t);

        if let Some((d, r)) = self_next.intersection(&other_next) {
            let n = d.norm();
            let elasticity = 0.975;
            let p = (1.0 + elasticity) * (self.v.dot(n) - other.v.dot(n)) / (self.m + other.m);
            let self_i = n * -p * self.m;
            let other_i = n * p * other.m;
            Some((
                self_i,
                other_i,
            ))
        } else {
            None
        }
    }
}

fn main() {
    let mut window = Window::new_flags(
        -1, -1,
        1024, 768,
        "Physics",
        &[WindowFlag::Async]
    ).unwrap();

    let window_w = window.width() as i32;
    let window_h = window.height() as i32;

    /*
    let G = 1.0E-2;
    let mut objects = Vec::with_capacity(64);
    for i in 0..objects.capacity() as i32 {
        let col = i % 8 - 4;
        let row = i / 8 - 4;
        objects.push(Object {
            p: Vec2::new(
                ((col * 64) + window_w / 2) as f64,
                ((row * 64) + window_h / 2) as f64,
            ),
            v: Vec2::new(
                0.0 * -(row as f64) / 96.0,
                0.0 * (col as f64) / 96.0,
            ),
            m: 1.0,
            r: 4.0,
        });
    }
    */

    let G = 6.6740831E-11; // G in m^3/(kg*s^2)
    let AU = 149_597_870_700.0; // AU in m
    let t = 3600.0; // hour in seconds
    let scale = (window.height() as f64) / (3.0 * AU);

    let mut objects = Vec::new();

    // Sun
    objects.push(Object {
        p: Vec2::new(
            0.0,
            0.0,
        ),
        v: Vec2::new(
            0.0,
            0.0,
        ),
        r: 6.9551E8, // radius in m
        m: 1.989E30, // mass in kg
        c: Color::rgb(0xFC, 0xD4, 0x40),
    });

    // Mercury
    objects.push(Object {
        p: Vec2::new(
            0.0,
            0.387 * AU,
        ),
        v: Vec2::new(
            47_362.0, // speed in m/s,
            0.0,
        ),
        r: 2.4397E6, // radius in m
        m: 3.3011E23, // mass in kg
        c: Color::rgb(0xA0, 0xA0, 0xA0),
    });

    // Venus
    objects.push(Object {
        p: Vec2::new(
            0.0,
            0.723 * AU,
        ),
        v: Vec2::new(
            36_020.0, // speed in m/s,
            0.0,
        ),
        r: 6.0518E6, // radius in m
        m: 4.867E24, // mass in kg
        c: Color::rgb(0xEE, 0xCB, 0x8B),
    });

    // Earth
    objects.push(Object {
        p: Vec2::new(
            0.0,
            AU,
        ),
        v: Vec2::new(
            29_800.0, // speed in m/s,
            0.0,
        ),
        r: 6.3781E6, // radius in m
        m: 5.972E24, // mass in kg
        c: Color::rgb(0x00, 0x77, 0xBE),
    });

    // Moon
    objects.push(Object {
        p: Vec2::new(
            0.0,
            AU + 3.844E8,
        ),
        v: Vec2::new(
            29_800.0 + 1023.0, // speed in m/s,
            0.0,
        ),
        r: 1.7371E6, // radius in m
        m: 7.34767309E22, // mass in kg
        c: Color::rgb(0xA0, 0xA0, 0xA0),
    });

    // Mars
    objects.push(Object {
        p: Vec2::new(
            0.0,
            1.523 * AU,
        ),
        v: Vec2::new(
            24_007.0, // speed in m/s,
            0.0,
        ),
        r: 3.3895E6, // radius in m
        m: 6.4171E23, // mass in kg
        c: Color::rgb(0xFF, 0x00, 0x00),
    });

    loop {
        window.clear();

        // a = F / m
        // F = G * m1 * m2 / (r * r)
        // a = G * m2 / (r * r)
        // First, calculate new velocities
        for i in 0..objects.len() {
            let p = {
                let object = &objects[i];
                object.p
            };

            let mut a = Vec2::new(0.0, 0.0);

            for j in 0..objects.len() {
                if j != i {
                    let other = &objects[j];
                    let d = other.p - p;
                    let r_squared = d.length_squared();
                    let a_mag = G * other.m / r_squared;
                    let theta = d.angle();
                    let (sin, cos) = theta.sin_cos();
                    a.x += a_mag * cos;
                    a.y += a_mag * sin;
                }
            }

            {
                let mut object = &mut objects[i];
                object.v += a * t;
            }
        }

        // Next, calculate collisions
        for i in 0..objects.len() {
            for j in 0..objects.len() {
                if j != i {
                    if let Some(impulse) = objects[i].collide(&objects[j], t) {
                        //TODO: Apply later to allow objects to collide at the same time?
                        objects[i].v += impulse.0;
                        objects[j].v += impulse.1;
                    }
                }
            }
        }


        for object in objects.iter() {
            let (current, left, right, next) = object.collision_geometry(t);
            //current.draw(&mut window, Color::rgb(0xFF, 0, 0), scale);
            next.draw(&mut window, object.c, scale);
        }

        // Finally, perform movements
        for object in objects.iter_mut() {
            object.p += object.v * t;
        }

        window.sync();

        for event in window.events() {
            match event.to_option() {
                EventOption::Quit(_) => return,
                _ => ()
            }
        }
    }
}
